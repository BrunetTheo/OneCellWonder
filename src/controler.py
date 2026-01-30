import cellStatus
import numpy as np
import interface

class Controller:
    def __init__(self,x:int,y:int,configFile:str,rulesFile:str):
        self.shape = (x,y)
        self.configFile = configFile
        self.rulesFile = rulesFile

        self.cellGrid = cellStatus.initialise_grid(self.rulesFile,self.configFile,x,y)
        self.interfce = interface.Interface((1000,1000),self)
    
    def update(self):
        self.cellGrid.update_grid()
    
    def getGrid(self):
        return self.cellGrid.getCellStatus().astype(int)

if __name__ == "__main__":
    c = Controller(100,100,"../confs/periodic/exempleCellConfig.txt","../confs/periodic/rules.txt")