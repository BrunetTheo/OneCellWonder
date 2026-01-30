import numpy as np
from dataclasses import dataclass
from src.parse_rules import AndRule
from scipy.signal import convolve2d
from parse_rules import *
from parse_cells import *
import utils as utils
import copy

@dataclass
class Cell:
    x: int
    y: int
    active_genes: np.array
    

def initialise_grid(name_file_rules,name_file_cell,X=20,Y=50,G=100):
    genes_rules,alive_rules = read_rules_file(name_file_rules)
    initial_cells = parse_cell_conf(name_file_cell)
    return CellGrid(X,Y,G,genes_rules=genes_rules,alive_rules=alive_rules,initial_cells=initial_cells)


class CellGrid:

    def __init__(self, X, Y, G, genes_rules,alive_rules,initial_cells=None, gene_names=None):
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
        self.genes_rules = genes_rules
        self.alive_rules = alive_rules

        # 1️⃣ Cellules vivantes / mortes
        self.cell_status = np.zeros((X, Y), dtype=bool)

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
            for gene in genes:
                    self.gene_content[x, y, gene] = 1

    def getCellStatus(self):
        return self.cell_status
    
    def _in_bounds(self, x, y):
        return 0 <= x < self.X and 0 <= y < self.Y


    def get_coords(self):
        """Return coordinates of all cells on the grid. Implies whether a cell is dead or alive"""
        pass
    
    def get_neighbors(self,n=1):
        """
        Docstring for get_neighbors
        
        n (int): number of wanted  
        """
        maskEven = utils.makeMask(True,n)        
        maskOdd = utils.makeMask(False,n)        
        
        to_int = np.array(self.cell_status,dtype=int)
        out_even = convolve2d(to_int, maskEven, mode="same")
        out_odd  = convolve2d(to_int, maskOdd,  mode="same")

        cols = np.arange(self.cell_status.shape[0])[None, :]
        evenCols = (cols % 2 == 0)

        neighbors = np.where(evenCols, out_even, out_odd) 
        return neighbors
        
    def neigboor_mask(self,applicable,n):
        maskEven = utils.makeMask(True,n)        
        maskOdd = utils.makeMask(False,n)        
        
        
        out_even = convolve2d(applicable, maskEven, mode="same")
        out_odd  = convolve2d(applicable, maskOdd,  mode="same")

        cols = np.arange(self.cell_status.shape[0])[None, :]
        evenCols = (cols % 2 == 0)

        neighbors = np.where(evenCols, out_even, out_odd) 
        return np.array(neighbors,dtype=bool)

        
    
    def apply_neighborhoodmask(self, gene_grid, rule_applied_grid, affected_neighborhood, gene_idx):
        """based on the affected_neighborhood, it creates a mask with radius affected_neighborhood, 
        and applies the mask to the origin of signal (found in rule_applied_grid)."""
        rows, cols = np.where(rule_applied_grid == 1)
        #get mask
        for r, c in zip(rows, cols):
            #check if y-coordinate is even or odd
            if c % 2 == 0:
                mask = utils.makeMask(True, affected_neighborhood)
            else:
                mask = utils.makeMask(False, affected_neighborhood)
            # apply mask
            gene_grid[r, c, gene_idx][mask == 1] = 1
        return gene_grid

    

    
    def get_genes(self, cell_indices=None):
        """Returns set of active genes per cell"""
        pass
    
    def validate_rule(self,positive_genes,negative_genes,neighboor_grid,n_neighboor=None):
        """ Check If a rule return true
        
        positive_genes = numpy array or gene that must be present e.g [1] when the second gene need to be there
        negatvie_genes = numpy array or gene that must be absent e.g [2,3] when the third and fourth gene must be absent

        """

        n_positive = len(positive_genes)
        if n_positive != 0:
            positive_gene_validation = np.sum(self.gene_content[:,:,positive_genes],axis=-1) == n_positive
        else:
            positive_gene_validation = np.ones_like(self.gene_content[:,:,0],dtype=bool)

        n_negative = len(negative_genes)
        if n_negative != 0:

            negative_gene_validation = np.sum(1-self.gene_content[:,:,negative_genes],axis=-1) == n_negative
        else:
            negative_gene_validation = np.ones_like(self.gene_content[:,:,0],dtype=bool)


        gene_validation = positive_gene_validation * negative_gene_validation

        if n_neighboor != None:
            gene_validation = gene_validation * (n_neighboor == neighboor_grid)
        else:
            gene_validation =  gene_validation * self.cell_status
        
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
        

    def propagate_genes(self):

        neighboor_grid = self.get_neighbors()


        new_genes = copy.deepcopy(self.gene_content)
        for rule in self.genes_rules:
            applicable = self.validate_rule(rule.positive_genes,rule.negative_genes,
                                           neighboor_grid=neighboor_grid,
                                           n_neighboor = rule.n_neighboor)
             
             
            extent = self.neigboor_mask(np.array(applicable,dtype=int),rule.propagation)

            new_genes[:,:,rule.active_gene] = new_genes[:,:,rule.active_gene] | extent

        self.gene_content = new_genes

    def create_alive_cell(self):
        neighboor_grid = self.get_neighbors()

        # First find potential new cells
        if len(self.alive_rules) == 0:
            return

        if len(self.alive_rules) >= 1:
            rule = self.alive_rules[0]
            new_alive = self.validate_rule(rule.positive_genes,rule.negative_genes,
                                           neighboor_grid=neighboor_grid,
                                           n_neighboor = rule.n_neighboor)
            if len(self.alive_rules ) > 1:
                rule = self.alive_rules[1]

                new_alive  = new_alive | self.validate_rule(rule.positive_genes,rule.negative_genes,
                                                            neighboor_grid=neighboor_grid,
                                                            n_neighboor = rule.n_neighboor)
        
        # remove allready alive

        new_alive = new_alive * (~self.cell_status)

        self.cell_status = self.cell_status | new_alive

    def update_grid(self):

        #compute N neighboor once:
        self.create_alive_cell()
        self.propagate_genes()
       
            #where the rules apply


        



