import numpy as np
from dataclasses import dataclass
from scipy.signal import convolve2d
from src.parse_rules import AndRule
from src.parse_rules import *
from src.parse_cells import *
import src.utils as utils
import copy


@dataclass
class Cell:
    x: int
    y: int
    active_genes: np.array


def initialise_grid(name_file_rules, name_file_cell, X=20, Y=50, G=100):
    genes_rules, alive_rules, ngene = read_rules_file(name_file_rules)
    initial_cells = parse_cell_conf(name_file_cell)
    return CellGrid(X, Y, ngene,
                    genes_rules=genes_rules,
                    alive_rules=alive_rules,
                    initial_cells=initial_cells)


class CellGrid:

    def __init__(self, X, Y, G, genes_rules, alive_rules,
                 initial_cells=None, gene_names=None):
        self.X = X
        self.Y = Y
        self.G = G
        self.genes_rules = genes_rules
        self.alive_rules = alive_rules

        self.cell_status = np.zeros((X, Y), dtype=int)
        self.gene_content = np.zeros((X, Y, G), dtype=int)

        if gene_names is None:
            self.gene_names = [f"gene_{i}" for i in range(G)]
        else:
            if len(gene_names) != G:
                raise ValueError("gene_names must have length G")
            self.gene_names = gene_names

        if initial_cells is not None:
            for cell in initial_cells:
                x, y = cell.x, cell.y
                if self._in_bounds(x, y):
                    self.cell_status[x, y] = 1
                for gene in cell.active_genes:
                    self.gene_content[x, y, gene] = 1

        self.propagate_genes()

    def getCellStatus(self):
        return self.cell_status

    def _in_bounds(self, x, y):
        return 0 <= x < self.X and 0 <= y < self.Y

    # ------------------------------------------------------------------
    # Convolution — all delegated to utils
    # ------------------------------------------------------------------

    def get_neighbors(self, n=1):
        """
        Count alive neighbours for each cell within radius n.
        Center cell is excluded from its own neighbour count.
        """
        return utils.adaptive_convolution(
            self.cell_status, n, self.X, self.Y, include_center=False
        )

    def neigboor_mask(self, applicable, n, candidate_coords=None):
        """
        Boolean mask: True where at least one neighbour within n is nonzero.

        candidate_coords : np.ndarray shape (N, 2), optional
            Superset of active coordinates in `applicable`. Passed through to
            sparse_convolution to avoid a redundant argwhere scan.
        """
        neighbors = utils.adaptive_convolution(
            applicable.copy(), n, self.X, self.Y,
            include_center=True,
            candidate_coords=candidate_coords
        )
        return neighbors.astype(bool)

    def inclusive_neigboor_mask(self, applicable, n, candidate_coords=None):
        """Boolean mask including the source cells themselves."""
        return (self.neigboor_mask(applicable, n, candidate_coords=candidate_coords)
                | applicable.astype(bool))

    # ------------------------------------------------------------------
    # Rule validation
    # ------------------------------------------------------------------

    def validate_rule(self, positive_genes, negative_genes,
                      neighboor_grid, potential_cell, n_neighboor=None):
        X, Y = self.gene_content.shape[:2]

        if len(positive_genes) > 0:
            positive_gene_validation = np.all(
                self.gene_content[:, :, positive_genes], axis=-1
            )
        else:
            positive_gene_validation = np.ones((X, Y), dtype=bool)

        if len(negative_genes) > 0:
            has_negative = np.any(
                self.gene_content[:, :, negative_genes], axis=-1
            )
            negative_gene_validation = ~has_negative
        else:
            negative_gene_validation = np.ones((X, Y), dtype=bool)

        gene_validation = positive_gene_validation & negative_gene_validation

        if n_neighboor is not None:
            gene_validation = gene_validation & (n_neighboor == neighboor_grid)
        else:
            gene_validation = gene_validation & potential_cell

        return gene_validation

    # ------------------------------------------------------------------
    # Gene propagation
    # ------------------------------------------------------------------

    def propagate_genes(self):
        # Mask neighbour count by cell_status: dead cells report 0 neighbours.
        # Without this, a dead cell with the right number of alive neighbours
        # could incorrectly match a rule with n_neighboor specified.
        neighboor_grid  = self.get_neighbors() * self.cell_status
        new_genes       = np.zeros_like(self.gene_content)
        expressed_genes = np.any(self.gene_content, axis=(0, 1))

        # Computed ONCE — applicable is always a subset of cell_status,
        # so alive_coords is a valid superset for all inclusive_neigboor_mask calls below.
        alive_coords = np.argwhere(self.cell_status != 0)

        for rule in self.genes_rules:
            if len(rule.positive_genes) > 0 and not np.all(expressed_genes[rule.positive_genes]):
                continue

            # For gene rules, only alive cells can express genes.
            # potential_cell (alive + neighbours) is only correct for create_alive_cell.
            applicable = self.validate_rule(
                rule.positive_genes, rule.negative_genes,
                neighboor_grid=neighboor_grid,
                potential_cell=self.cell_status,
                n_neighboor=rule.n_neighboor
            )

            extent = self.inclusive_neigboor_mask(
                np.array(applicable, dtype=int), rule.propagation,
                candidate_coords=alive_coords  # avoids argwhere inside sparse_convolution
            )

            new_genes[:, :, rule.active_gene] = (
                new_genes[:, :, rule.active_gene] | extent
            )

        self.gene_content = new_genes

    # ------------------------------------------------------------------
    # Cell creation
    # ------------------------------------------------------------------

    def create_alive_cell(self):
        if len(self.alive_rules) == 0:
            return

        neighboor_grid = self.get_neighbors()
        potential_cell = self.inclusive_neigboor_mask(self.cell_status.copy(), 1)

        rule = self.alive_rules[0]
        new_alive = self.validate_rule(
            rule.positive_genes, rule.negative_genes,
            neighboor_grid=neighboor_grid,
            potential_cell=potential_cell,
            n_neighboor=rule.n_neighboor
        )

        for rule in self.alive_rules[1:]:
            new_alive |= self.validate_rule(
                rule.positive_genes, rule.negative_genes,
                neighboor_grid=neighboor_grid,
                potential_cell=potential_cell,
                n_neighboor=rule.n_neighboor
            )

        new_alive = new_alive * (1 - self.cell_status)
        self.cell_status = self.cell_status | new_alive

    def update_grid(self):
        self.create_alive_cell()
        self.propagate_genes()