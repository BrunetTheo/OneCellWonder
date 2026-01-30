import numpy as np
from dataclasses import dataclass


@dataclass
class AndRule:
    positive_genes : np.array
    negative_genes  : np.array
    n_neighboor : int
    propagation : int
    

class CellGrid:
    def __init__(self, state_grid, gene_grid):
        """
        C: Number of Cells
        G: Number of genes
        alive_grid: np.array [grid dimension]
        gene_grid: np.array [W x L x G]
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
    
    def validate_rule(self,positive_genes,negative_genes,n_neighboor=None):
        """ Check If a rule return true
        
        positive_genes = numpy array or gene that must be present e.g [1] when the second gene need to be there
        negatvie_genes = numpy array or gene that must be absent e.g [2,3] when the third and fourth gene must be absent

        """

        n_positive = len(positive_genes)
        if n_positive != 0:
            positive_gene_validation = np.sum(self.gene_grid[:,:,positive_genes],axis=-1) == n_positive
        else:
            positive_gene_validation = np.ones_like(self.gene_grid[:,:,0],dtype=bool)

        n_negative = len(negative_genes)
        if n_negative != 0:

            negative_gene_validation = np.sum(1-self.gene_grid[:,:,negative_genes],axis=-1) == n_negative
        else:
            negative_gene_validation = np.ones_like(self.gene_grid[:,:,0],dtype=bool)


        gene_validation = positive_gene_validation * negative_gene_validation

        if n_neighboor != None:
            gene_validation = gene_validation * (n_neighboor == self.get_n_neighboor())
        
        return gene_validation


    def match_rules(self, rules):
        """matches the rules to cell status """
        if rules != []:
            rule = rules[0]
            init = self.validate_rule(rule.positive_genes,rule.negative_genes,rule.n_neighboor)
            for rule in rules[1:]:
                init *= self.validate_rule(rule.positive_genes,rule.negative_genes,rule.n_neighboor)
            return  np.zeros_like(self.gene_grid[:,:,0],dtype=bool)
        else:
            return 
    def update_grid(self):
        pass



