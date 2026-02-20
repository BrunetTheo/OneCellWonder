import numpy as np
from dataclasses import dataclass
from src.parse_rules import AndRule
from src.parse_rules import *
from src.parse_cells import *
import src.utils as utils


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

        self.cell_status = np.zeros((X, Y), dtype=np.int8)

        # gene_content: (G, X, Y) int8 array.
        # Gene g is at gene_content[g] — a 2D (X, Y) slice.
        self.gene_content = np.zeros((G, X, Y), dtype=np.int8)

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
                    self.gene_content[gene, x, y] = 1

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
        """Boolean mask: True where at least one neighbour within n is nonzero."""
        neighbors = utils.adaptive_convolution(
            applicable, n, self.X, self.Y,
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
        """
        Returns a (X, Y) boolean mask of cells where the rule applies.
        gene_content[g] is a contiguous (X, Y) slice — same cost as dict lookup.
        """
        if len(positive_genes) > 0:
            validation = self.gene_content[positive_genes[0]].astype(bool)
            for g in positive_genes[1:]:
                validation &= self.gene_content[g].astype(bool)
        else:
            validation = np.ones((self.X, self.Y), dtype=bool)

        for g in negative_genes:
            validation &= ~self.gene_content[g].astype(bool)

        if n_neighboor is not None:
            validation &= (n_neighboor == neighboor_grid)
        else:
            validation &= potential_cell.astype(bool)

        return validation

    # ------------------------------------------------------------------
    # Gene propagation
    # ------------------------------------------------------------------

    def propagate_genes(self):
        # Mask neighbour count by cell_status: dead cells report 0 neighbours.
        neighboor_grid = self.get_neighbors() * self.cell_status

        # expressed_genes: one vectorised reduction over (G, X, Y).
        expressed_genes = set(np.where(self.gene_content.any(axis=(1, 2)))[0])

        # Computed ONCE — applicable is always a subset of cell_status.
        alive_coords = np.argwhere(self.cell_status != 0)

        new_genes = np.zeros((self.G, self.X, self.Y), dtype=np.int8)

        for rule in self.genes_rules:
            if len(rule.positive_genes) > 0 and not all(
                g in expressed_genes for g in rule.positive_genes
            ):
                continue

            applicable = self.validate_rule(
                rule.positive_genes, rule.negative_genes,
                neighboor_grid=neighboor_grid,
                potential_cell=self.cell_status,
                n_neighboor=rule.n_neighboor
            )

            extent = self.inclusive_neigboor_mask(
                applicable.astype(np.int8), rule.propagation,
                candidate_coords=alive_coords
            )

            new_genes[rule.active_gene] |= extent

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

        new_alive = new_alive & ~self.cell_status.astype(bool)
        self.cell_status = (self.cell_status | new_alive).astype(np.int8)

    def update_grid(self):
        self.create_alive_cell()
        self.propagate_genes()