import re
import numpy as np
from dataclasses import dataclass


@dataclass
class AndRule:
    positive_genes : np.array
    negative_genes  : np.array
    n_neighboor : int
    propagation : int
    active_gene: int

def parse_and_rule(rule_str: str, active_gene: int) -> AndRule:
    """
    Parse a single rule like:
    [1,2,not(3),n(6)]7
    """
    # Extract content inside brackets and propagation
    match = re.fullmatch(r"\[(.*?)\](\d+)", rule_str.strip())
    if not match:
        raise ValueError(f"Invalid rule format: {rule_str}")

    content, propagation = match.groups()
    propagation = int(propagation)

    positive = []
    negative = []
    n_neighboor = None

    for token in content.split(","):
        token = token.strip()

        if token.startswith("not("):
            negative.append(int(token[4:-1]))
        elif token.startswith("n("):
            n_neighboor = int(token[2:-1])
        else:
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

def read_rules_file(filepath: str) -> List[AndRule]:
    """
    Read a rules file and return a list of AndRule objects.
    Each line corresponds to one active_gene (line index).
    """
    all_rules: List[AndRule] = []

    with open(filepath, "r") as f:
        for active_gene, line in enumerate(f):
            line = line.strip()

            # Skip empty lines and comments
            if not line or line.startswith("#"):
                continue

            rules = parse_rule_line(line, active_gene)
            all_rules.extend(rules)

    return all_rules


if __name__ == "__main__":
    line = "[1,2,not(3),n(6),4,not(1)]7 || [1,n(2)]5"
    rules = parse_rule_line(line, active_gene=1)

    for r in rules:
        print(r)
