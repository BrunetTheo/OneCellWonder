import numpy as np
from dataclasses import dataclass
from src.parse_rules import AndRule
from scipy.signal import convolve2d
from src.parse_rules import *
from src.parse_cells import *
import src.utils as utils
import copy

@dataclass
class Cell:
    x: int
    y: int
    active_genes: np.array
    

def initialise_grid(name_file_rules,name_file_cell,X=20,Y=50,G=100):
    genes_rules,alive_rules,ngene = read_rules_file(name_file_rules)
    initial_cells = parse_cell_conf(name_file_cell)
    return CellGrid(X,Y,ngene,genes_rules=genes_rules,alive_rules=alive_rules,initial_cells=initial_cells)


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
                for gene in genes:
                        self.gene_content[x, y, gene] = 1
        
        self.propagate_genes()

    def getCellStatus(self):
        return self.cell_status
    
    def _in_bounds(self, x, y):
        return 0 <= x < self.X and 0 <= y < self.Y


    
    def get_neighbors(self,n=1):
        """
        Docstring for get_neighbors
        
        n (int): number of wanted  
        """
        maskEven = utils.makeMask(True,n)        
        maskOdd = utils.makeMask(False,n)        
        
        out_even = convolve2d(self.cell_status, maskEven, mode="same")
        out_odd  = convolve2d(self.cell_status, maskOdd,  mode="same")

        cols = np.arange(self.cell_status.shape[1])[None, :]
        evenCols = (cols % 2 == 0)

        neighbors = np.where(evenCols, out_even, out_odd) 
        return neighbors


    def sparse_convolution(self, matrix, n):
        """Only convolve around non-zero cells - CORRECTED."""
        
        active_coords = np.argwhere(matrix != 0)
        
        if len(active_coords) == 0:
            return np.zeros_like(matrix)
        
        result = np.zeros_like(matrix, dtype=int)

        iseven_for_even_cols = True
        iseven_for_odd_cols = False
        
        if n % 2 != 0:
            iseven_for_even_cols, iseven_for_odd_cols = iseven_for_odd_cols, iseven_for_even_cols
        
        # Get masks with swapped parameters
        maskEven_with_center = utils.makeMask(iseven_for_even_cols, n, include_center=True)
        maskOdd_with_center = utils.makeMask(iseven_for_odd_cols, n, include_center=True)
        
        
        for x, y in active_coords:
            mask = maskEven_with_center if y % 2 == 0 else maskOdd_with_center
            
            x_start = max(0, x - n)
            x_end = min(self.X, x + n + 1)
            y_start = max(0, y - n)
            y_end = min(self.Y, y + n + 1)
            
            mask_x_start = n - (x - x_start)
            mask_x_end = n + (x_end - x)
            mask_y_start = n - (y - y_start)
            mask_y_end = n + (y_end - y)
            
            mask_slice = mask[mask_x_start:mask_x_end, mask_y_start:mask_y_end]
            
            # CORRECTED: Use matrix[x, y] not neighborhood
            result[x_start:x_end, y_start:y_end] += matrix[x, y] * mask_slice
        
        return result
    

    def convolution(self, matrix, n):
      
        iseven_for_even_cols = True
        iseven_for_odd_cols = False

        if n % 2 != 0:
            iseven_for_even_cols, iseven_for_odd_cols = iseven_for_odd_cols, iseven_for_even_cols
        
        # Get masks with swapped parameters
        maskEven = utils.makeMask(iseven_for_even_cols, n, include_center=True)
        maskOdd = utils.makeMask(iseven_for_odd_cols, n, include_center=True)

        # Create proper copies and ensure we don't modify the input
        even = matrix.copy()  # Use .copy() for numpy arrays
        odd = matrix.copy()
        
        even[:, 1::2] = 0
        odd[:, ::2] = 0
        
        out_even = convolve2d(even, maskEven, mode="same")
        out_odd = convolve2d(odd, maskOdd, mode="same")    
        
        neighbors = out_even + out_odd
        return neighbors
    
    def adaptive_convolution(self, matrix, n):
        """Choose sparse or dense convolution based on occupancy."""
        occupancy = np.count_nonzero(matrix) / matrix.size
        
        if occupancy < 0.1:  # Less than 10% occupied
            return self.sparse_convolution(matrix, n)
        else:
            return self.convolution(matrix, n)  # Original dense version

    def neigboor_mask(self, applicable, n):
        # Create a true independent copy before passing to convolution
        neighbors = self.adaptive_convolution(applicable.copy(), n)
        return np.array(neighbors, dtype=bool)
    

    def inclusive_neigboor_mask(self,applicable,n):
        return self.neigboor_mask(applicable,n) | applicable
        
        
    
    def validate_rule(self, positive_genes, negative_genes, 
                 neighboor_grid, potential_cell, n_neighboor=None):
        """ Check If a rule return true
        
        positive_genes = numpy array or gene that must be present e.g [1] when the second gene need to be there
        negatvie_genes = numpy array or gene that must be absent e.g [2,3] when the third and fourth gene must be absent

        """
        X, Y = self.gene_content.shape[:2]
        
        # Positive genes: use np.all (faster)
        if len(positive_genes) > 0:
            positive_gene_validation = np.all(
                self.gene_content[:, :, positive_genes], axis=-1
            )
        else:
            positive_gene_validation = np.ones((X, Y), dtype=bool)
        
        # Negative genes: check if ANY are present, then invert
        if len(negative_genes) > 0:
            has_negative = np.any(
                self.gene_content[:, :, negative_genes], axis=-1
            )
            negative_gene_validation = ~has_negative
        else:
            negative_gene_validation = np.ones((X, Y), dtype=bool)
        
        # Use bitwise AND (faster than multiplication)
        gene_validation = positive_gene_validation & negative_gene_validation
        
        if n_neighboor is not None:
            gene_validation = gene_validation & (n_neighboor == neighboor_grid)
        else:
            gene_validation = gene_validation & potential_cell
        
        return gene_validation


    def propagate_genes(self):

        neighboor_grid = self.get_neighbors()
        potential_cell = self.inclusive_neigboor_mask(self.cell_status.copy(),1) #Cells and their one radius neigboor
        new_genes = np.zeros_like(self.gene_content)

        expressed_genes = np.any(self.gene_content, axis=(0, 1))  # shape (G,) — one bool per gene

        for rule in self.genes_rules:
            # If any required positive gene is absent from the entire grid,
            # validate_rule can never return True anywhere — skip it entirely
            if len(rule.positive_genes) > 0 and not np.all(expressed_genes[rule.positive_genes]):
                continue
            applicable = self.validate_rule(rule.positive_genes,rule.negative_genes,
                                           neighboor_grid=neighboor_grid,
                                           potential_cell=potential_cell,
                                           n_neighboor = rule.n_neighboor)
             
            applicable = applicable * self.cell_status   #Genes diffuses only from alive cell
            extent = self.inclusive_neigboor_mask(np.array(applicable,dtype=int),rule.propagation) 
            

            new_genes[:,:,rule.active_gene] = new_genes[:,:,rule.active_gene] | extent

        self.gene_content = new_genes

    def create_alive_cell(self):
        neighboor_grid = self.get_neighbors()
        potential_cell = self.inclusive_neigboor_mask(self.cell_status.copy(),1) #Cells and their one radius neigboor

        # First find potential new cells
        if len(self.alive_rules) == 0:
            return

        if len(self.alive_rules) >= 1:
            rule = self.alive_rules[0]
            new_alive = self.validate_rule(rule.positive_genes,rule.negative_genes,
                                           neighboor_grid=neighboor_grid,
                                            potential_cell=potential_cell,
                                           n_neighboor = rule.n_neighboor)
            if len(self.alive_rules ) > 1:
                for rule in  self.alive_rules[1:]:
                    new_alive  = new_alive | self.validate_rule(rule.positive_genes,rule.negative_genes,
                                                            neighboor_grid=neighboor_grid,
                                                            potential_cell=potential_cell,
                                                            n_neighboor = rule.n_neighboor)
        
        # remove allready alive

        new_alive = new_alive * (1 - self.cell_status)
        #return new_alive
        self.cell_status = self.cell_status | new_alive

    def update_grid(self):

        #compute N neighboor once:
      
        self.create_alive_cell()
        self.propagate_genes()
       
            #where the rules apply