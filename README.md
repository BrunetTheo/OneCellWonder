# OneCellWonder



# Install
mamba create -n  one_cell_wonder conda-forge::pygame numpy scipy



# Run
```
python main.py confs/default/ --X 80 --Y 40
```



# System and Engine Overview

This is a gene-based cellular automata system where cells exist on a 2D hexagonal grid and contain genetic information that determines their behavior and propagation.

## Initialization

The system is initialized with the following parameters:

- **Grid dimensions**: `X × Y` spatial grid (default: 20 × 50)
- **Genes**: `G` different genes per cell (default: 3)
- **Cell status**: Binary grid `cell_status[x, y] ∈ {0, 1}` indicating alive/dead cells
- **Gene content**: 3D array `gene_content[x, y, g] ∈ {0, 1}` indicating which genes are active in each cell

**Initial conditions:**
1. Place initial living cells at specified coordinates
2. Assign active genes to these initial cells
3. Perform initial gene propagation to nearby cells

## Update Equations

The system evolves through two main operations executed each timestep:

### 1. Cell Birth (`create_alive_cell`)

New cells are born based on rules that check:
- **Gene conditions**: Presence of required genes (positive) and absence of forbidden genes (negative)
- **Neighbor count**: Number of living neighbors `N[x, y]`

For each alive rule with conditions (P⁺, P⁻, n):

```
valid[x, y] = (Σ gene_content[x, y, g] == |P⁺| for g ∈ P⁺) 
            ∧ (Σ (1 - gene_content[x, y, g]) == |P⁻| for g ∈ P⁻)
            ∧ (N[x, y] == n)
```

New cells are born where any rule is satisfied and no cell currently exists:

```
new_alive = (⋁ valid_rule[x, y]) ∧ ¬cell_status[x, y]
cell_status[x, y] ← cell_status[x, y] ∨ new_alive[x, y]
```

### 2. Gene Propagation (`propagate_genes`)

Genes spread to neighboring cells based on propagation rules:

For each gene propagation rule (P⁺, P⁻, n, g_target, r):

1. **Find source cells** that satisfy the rule:
```
applicable[x, y] = (gene conditions satisfied) ∧ (N[x, y] == n)
```

2. **Create neighborhood mask** of radius `r` around each applicable cell

3. **Activate target gene** in all cells within the mask:
```
gene_content[x, y, g_target] ← gene_content[x, y, g_target] ∨ mask[x, y]
```

### Hexagonal Neighbor Counting

The system uses hexagonal grid topology where neighbor counts depend on column parity:
- **Even columns**: Use even-column hexagonal mask
- **Odd columns**: Use odd-column hexagonal mask

```
N[x, y] = Σ cell_status[neighbors(x, y)]
```

where `neighbors(x, y)` respects hexagonal adjacency rules.

## Key Features

- Gene-based rules control both cell birth and gene propagation
- Hexagonal grid topology for realistic spatial interactions
- Genes can spread beyond living cells to influence future cell birth
- Multiple rules can be combined with AND/OR logic