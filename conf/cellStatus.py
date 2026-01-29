import numpy as np

class CellGrid:
    def __init__(self, coords, neighbors, genes, gene_names):
        """
        N: Number of cells on grid 
        G: Number of genes
        coords: np.array [N,2] of x,y positions
        neighbors: np.array [N] number of neighbors per cell
        genes: np.array [N,G] boolean matrix for gene presence
        gene_names: list of gene names corresponding to columns of genes
        """
        self.coords = coords
        self.neighbors = neighbors
        self.genes = genes
        self.gene_names = gene_names
        # Map gene name -> column index for convenience
        self.gene_idx = {name: i for i, name in enumerate(gene_names)}

    def get_coords(self, grid):
        """Return coordinates of all cells or selected cells. Implies whether a cell is dead or alive"""
        
        pass
       
    def get_genes(self, cell_indices=None):
        """Returns set of active genes per cell"""
        pass
       

    def match_rule(self, rule):
        """matches the rules to cell status """
        pass
