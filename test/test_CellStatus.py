import unittest
import cellStatus
import parse_rules 

class TestParseRules(unittest.TestCase):\
    def testParseAndRule(self):
        test_1 = "[0]"
        output = parse_rules.parse_and_rule(test_1):
        self.assertisinstance(output,parse_rules.AndRule)
        self.assertequal()

    def testParseRuleLine(self):

    def testReadRulesFile(self):
        pass
        









class TestCellGrid(unittest.TestCase):
    def setUp(self):
        CellGrid(X=5,
                 Y=5,
                 G=0,
                 genes_rules='',
                 alive_rules='',
                 initial_cells='')

    def test_init():
        
class CellGrid:

    def __init__(self, X, Y, G, genes_rules,alive_rules,initial_cells=None, gene_names=None):
