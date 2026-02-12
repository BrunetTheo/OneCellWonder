import numpy as np
import math
from scipy.signal import convolve2d

def makeMask(iseven,n):
    """
    Create a mask of size n
    
    iseven (bool): if the center of the mask is on an even x or not
    n (int):  of the neighbourhood 
    """
    meven=np.array([[1,1,1],[1,0,1],[0,1,0]])
    modd = np.array([[0,1,0],[1,0,1],[1,1,1]])
    
    mask = np.zeros((n*2+1,n*2+1))
    mask[n,n] = 1
    for i in range(n):

        center = meven if iseven else modd
        notCenter = modd if iseven else meven
        
        out_even = convolve2d(mask, center, mode="same")
        out_odd  = convolve2d(mask, notCenter,  mode="same")

        rows = np.arange(mask.shape[0])[None, :]
        mask_even = (rows % 2 == 0)

        mask = np.where(mask_even, out_even, out_odd) 
    mask = mask>0
    mask[n,n] = 0 #Remove the original cell
    return mask

