import cellStatus
import numpy as np
import interface

class Controller:
    def __init__(self,x:int,y:int,configFile:str,rulesFile:str):
        self.shape = (x,y)
        self.configFile = configFile
        self.rulesFile = rulesFile
        self.show = -1
        self.cellGrid = cellStatus.initialise_grid(self.rulesFile,self.configFile,x,y)
        self.interfce = interface.Interface((1000,1000),self)
    
    def update(self):
        self.cellGrid.update_grid()
    
    def getGrid(self):
        if self.show == -1:
            return self.cellGrid.getCellStatus().astype(int)
        if self.show == 0:
            return (self.cellGrid.gene_content[:,:,1]).astype(int)
        if self.show == 1:
            return (self.cellGrid.gene_content[:,:,2]).astype(int)

if __name__ == "__main__":
    c = Controller(100,100,"../confs/periodic/exempleCellConfig.txt","../confs/periodic/rules.txt")