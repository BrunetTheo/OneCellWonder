import numpy as np

class CellGrid:
    def __init__(self, state_grid, gene_grid):
        """
        C: Number of Cells
        G: Number of genes
        state_grid: np.array [grid dimension]
        gene_grid: np.array [G x C]
        """
        self.state_grid = state_grid
        self.gene_grid = gene_grid

    def get_coords(self):
        """Return coordinates of all cells on the grid. Implies whether a cell is dead or alive"""
        pass
    def get_neighbors(self):
        pass 
    
    
    def apply_neighborhoodmask(self, gene_grid, rule_applied_grid, affected_neighborhood, gene_idx):
        """based on the affected_neighborhood, it creates a mask of size affected_neighborhood, 
        and adds applies the mask to the origin of signal, which is indicated in rule_applied_grid rule_applied_grid: 
        matrix of coordinates that tells you where signal origin lies"""
        rows, cols = gene_grid
        mask = np.zeros(gene_grid, dtype= bool)

        # use a mask per cell of origin
            
            

    


    
    def get_genes(self, cell_indices=None):
        """Returns set of active genes per cell"""
        pass
       

    def match_rule(self, rule):
        """matches the rules to cell status """
        pass
    
    
    def update_grid(self):
        pass
    
    




# 