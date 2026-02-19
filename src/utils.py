import numpy as np
from scipy.signal import convolve2d
from functools import lru_cache

try:
    from numba import njit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    def njit(fn):
        return fn


# ---------------------------------------------------------------------------
# Mask construction
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def makeMask(iseven, n, include_center=False):
    meven = np.array([[1, 1, 1],[1, 0, 1],[0, 1, 0]])
    modd  = np.array([[0, 1, 0],[1, 0, 1],[1, 1, 1]])
    mask = np.zeros((n * 2 + 1, n * 2 + 1))
    mask[n, n] = 1
    for _ in range(n):
        center    = meven if iseven else modd
        notCenter = modd  if iseven else meven
        out_even  = convolve2d(mask, center,    mode="same")
        out_odd   = convolve2d(mask, notCenter, mode="same")
        rows      = np.arange(mask.shape[0])[None, :]
        mask      = np.where(rows % 2 == 0, out_even, out_odd)
    mask = mask > 0
    mask[n, n] = False
    if include_center:
        mask[n, n] = True
    return mask


# ---------------------------------------------------------------------------
# Dense convolution
# ---------------------------------------------------------------------------

def dense_convolution(matrix, n, X, Y, include_center=True):
    iseven   = n % 2 == 0
    maskEven = makeMask(iseven,     n, include_center=include_center)
    maskOdd  = makeMask(not iseven, n, include_center=include_center)
    even = matrix.copy(); even[:, 1::2] = 0
    odd  = matrix.copy(); odd[:, ::2]  = 0
    return (convolve2d(even, maskEven, mode="same") +
            convolve2d(odd,  maskOdd,  mode="same"))


# ---------------------------------------------------------------------------
# Sparse convolution (Numba kernel + Python fallback)
# ---------------------------------------------------------------------------

@njit
def _sparse_kernel(matrix, result, active_coords, mask_even, mask_odd, X, Y, n):
    for idx in range(len(active_coords)):
        x   = active_coords[idx, 0]
        y   = active_coords[idx, 1]
        val = matrix[x, y]
        x_start = x - n if x - n > 0 else 0
        x_end   = x + n + 1 if x + n + 1 < X else X
        y_start = y - n if y - n > 0 else 0
        y_end   = y + n + 1 if y + n + 1 < Y else Y
        mx_start = n - (x - x_start);  mx_end = n + (x_end - x)
        my_start = n - (y - y_start);  my_end = n + (y_end - y)
        mask = mask_even if y % 2 == 0 else mask_odd
        for i in range(x_end - x_start):
            for j in range(y_end - y_start):
                result[x_start + i, y_start + j] += val * mask[mx_start + i, my_start + j]
    return result


def sparse_convolution(matrix, n, X, Y, include_center=True, candidate_coords=None):
    """
    Sparse hex-grid convolution. Only processes non-zero cells.

    Parameters
    ----------
    candidate_coords : np.ndarray shape (N, 2), optional
        Pre-computed superset of active coordinates (e.g. alive_coords).
        Active coords are derived by filtering candidates where matrix != 0,
        avoiding a full np.argwhere scan over the whole grid.
        If None, falls back to np.argwhere.
    """
    if candidate_coords is not None and len(candidate_coords) > 0:
        active_mask   = matrix[candidate_coords[:, 0], candidate_coords[:, 1]] != 0
        active_coords = candidate_coords[active_mask]
    else:
        active_coords = np.argwhere(matrix != 0)

    if len(active_coords) == 0:
        return np.zeros_like(matrix, dtype=int)

    iseven    = n % 2 == 0
    mask_even = makeMask(iseven,     n, include_center=include_center).astype(np.int64)
    mask_odd  = makeMask(not iseven, n, include_center=include_center).astype(np.int64)

    if NUMBA_AVAILABLE:
        result = np.zeros((X, Y), dtype=np.int64)
        return _sparse_kernel(
            matrix.astype(np.int64), result,
            active_coords.astype(np.int64),
            mask_even, mask_odd,
            np.int64(X), np.int64(Y), np.int64(n),
        )
    else:
        result = np.zeros_like(matrix, dtype=int)
        for x, y in active_coords:
            mask = mask_even if y % 2 == 0 else mask_odd
            x_start = max(0, x - n);  x_end = min(X, x + n + 1)
            y_start = max(0, y - n);  y_end = min(Y, y + n + 1)
            mx_start = n - (x - x_start);  mx_end = n + (x_end - x)
            my_start = n - (y - y_start);  my_end = n + (y_end - y)
            result[x_start:x_end, y_start:y_end] += (
                matrix[x, y] * mask[mx_start:mx_end, my_start:my_end]
            )
        return result


# ---------------------------------------------------------------------------
# Adaptive convolution
# ---------------------------------------------------------------------------

def adaptive_convolution(matrix, n, X, Y, include_center=True, threshold=0.1,
                         candidate_coords=None):
    """
    Choose sparse or dense convolution based on matrix occupancy.

    candidate_coords is passed through to sparse_convolution to avoid
    redundant argwhere calls in the caller.
    """
    occupancy = np.count_nonzero(matrix) / matrix.size
    if occupancy < threshold:
        return sparse_convolution(matrix, n, X, Y,
                                  include_center=include_center,
                                  candidate_coords=candidate_coords)
    else:
        return dense_convolution(matrix, n, X, Y, include_center=include_center)