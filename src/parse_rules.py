import re
import numpy as np
from dataclasses import dataclass
from typing import List, Tuple


@dataclass
class AndRule:
    positive_genes : np.array
    negative_genes  : np.array
    n_neighboor : int
    propagation : int
    active_gene: int

def parse_and_rule(rule_str: str, active_gene: int) -> AndRule:
    """
    Supports:
      [1,2,not(3),n(6)]7
      [n(1)]            # propagation defaults to 0
    """
    match = re.fullmatch(r"\[(.*?)\](\d*)", rule_str.strip())
    if not match:
        raise ValueError(f"Invalid rule format: {rule_str}")

    content, propagation = match.groups()
    propagation = int(propagation) if propagation else 1

    positive = []
    negative = []
    n_neighboor = None

    if content.strip():
        for token in content.split(","):
            token = token.strip()

            if token.startswith("not(") and token.endswith(")"):
                negative.append(int(token[4:-1]))

            elif token.startswith("n(") and token.endswith(")"):
                n_neighboor = int(token[2:-1])

            else:
                if not token.isdigit():
                    raise ValueError(f"Invalid token '{token}' in rule '{rule_str}'")
                positive.append(int(token))

    return AndRule(
        positive_genes=np.array(positive, dtype=int),
        negative_genes=np.array(negative, dtype=int),
        n_neighboor=n_neighboor,
        propagation=propagation,
        active_gene=active_gene,
    )

def parse_rule_line(line: str, active_gene: int):
    """
    Parse a full line like:
    [1,2,not(3),n(6)]7 || [1]5
    """
    rule_strings = line.split("||")
    return [parse_and_rule(r, active_gene) for r in rule_strings]


from typing import List

def read_rules_file(filepath: str) -> Tuple[List[AndRule], List[AndRule]]:
    """
    Read a rules file.

    Returns:
        gene_rules:  AndRule from all lines except the last one
        alive_rules: AndRule from the last line to create new cell (active_gene = -1)
    """
    with open(filepath, "r") as f:
        lines = [
            line.strip()
            for line in f
            if line.strip() and not line.strip().startswith("#")
        ]

    if not lines:
        return [], []

    gene_rules: List[AndRule] = []
    alive_rules: List[AndRule] = []

    # All lines except the last → normal active_gene
    for active_gene, line in enumerate(lines[:-1]):
        gene_rules.extend(parse_rule_line(line, active_gene))

    # Last line → active_gene = -1
    alive_rules.extend(parse_rule_line(lines[-1], active_gene=-1))

    return gene_rules, alive_rules




if __name__ == "__main__":
    line = "[1,2,not(3),n(6),4,not(1)]7 || [1,n(2)]5"
    rules = parse_rule_line(line, active_gene=1)

    for r in rules:
        print(r)
