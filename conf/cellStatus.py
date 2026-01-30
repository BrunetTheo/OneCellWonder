import numpy as np
from dataclasses import dataclass


@dataclass
class AndRule:
    positive_genes : np.array
    negative_genes  : np.array
    n_neighboor : int
    propagation : int
    active_gene: int

@dataclass
class Cell:
    x: int
    y: int
    active_genes: np.array
    

class CellGrid:

    def __init__(self, X, Y, G, initial_cells=None, gene_names=None):
        """
        Grille spatiale de cellules avec contenu génétique.

        Paramètres
        ----------
        X : int
            Taille spatiale en x
        Y : int
            Taille spatiale en y
        G : int
            Nombre de gènes
        initial_cells : list of tuples (x,y), optionnel
            Positions initiales des cellules vivantes
        initial_genes : list of tuples (x,y,g), optionnel
            Gènes présents initialement dans les cellules
        gene_names : list of str, optionnel
            Noms des gènes
        """

        self.X = X
        self.Y = Y
        self.G = G

        # 1️⃣ Cellules vivantes / mortes
        self.cell_status = np.zeros((X, Y), dtype=int)

        # 2️⃣ Contenu génétique
        self.gene_content = np.zeros((X, Y, G), dtype=int)

        # Noms des gènes
        if gene_names is None:
            self.gene_names = [f"gene_{i}" for i in range(G)]
        else:
            if len(gene_names) != G:
                raise ValueError("gene_names doit être de longueur G")
            self.gene_names = gene_names

        # Initialisation des cellules vivantes
        if initial_cells is not None:
            for cell in initial_cells:
                x, y = cell.x,cell.y
                if self._in_bounds(x, y):
                    self.cell_status[x, y] = 1

                genes = cell.active_genes
                for gene in genes
                    self.gene_content[x, y, g] = 1


    def _in_bounds(self, x, y):
        return 0 <= x < self.X and 0 <= y < self.Y


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



