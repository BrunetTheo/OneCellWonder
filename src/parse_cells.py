from ast import literal_eval
from dataclasses import dataclass
import numpy as np

@dataclass
class Cell:
    x: int
    y: int
    active_genes: np.array

def parse_cell_conf(txt_path):
    cells_list = []
    with open(txt_path, 'r') as txt:
        for line in txt:
            elements = line.split(';')
            cells_list.append(Cell(
                x=int(elements[0]),
                y=int(elements[1]),
                active_genes=np.asarray(literal_eval(elements[2]))
            ))
    return cells_list
