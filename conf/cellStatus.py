import numpy as np

class CellGrid:
    def __init__(self, alive_grid, gene_grid):
        """
        C: Number of Cells
        G: Number of genes
        alive_grid: np.array [grid dimension]
        gene_grid: np.array [G x C]
        """
        self.alive_grid = alive_grid
        self.gene_grid = gene_grid

    def get_coords(self):
        """Return coordinates of all cells on the grid. Implies whether a cell is dead or alive"""
        pass
       
    def get_genes(self, cell_indices=None):
        """Returns set of active genes per cell"""
        pass
       

    def match_rule(self, rule):
        """matches the rules to cell status """
        pass
    def update_grid(self):
        pass
