import math
from HiggsAnalysis.CombinedLimit.PhysicsModel import *

### This is the base python class to study the SpinZero structure

class HighmassModel(PhysicsModel):
    def __init__(self):
        self.mHRange = []
        self.muAsPOI = False
        self.muFloating = False
        self.meanAsPOI = False
        self.widthAsPOI = False
        self.pois = {}
        self.verbose = False
    def setModelBuilder(self, modelBuilder):
        PhysicsModel.setModelBuilder(self,modelBuilder)
        self.modelBuilder.doModelBOnly = False

    def getYieldScale(self,bin,process):
        "Split in production and decay, and call getHiggsSignalYieldScale; return 1 for backgrounds "
        self.my_norm = 1 
        print "Process {0} will scale by {1}".format(process,self.my_norm)
        return self.my_norm
            

    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if 'muAsPOI' in po: 
                print "Will consider the signal strength as a parameter of interest"
                self.muAsPOI = True
                self.muFloating = True
            if 'mwAsPOI' in po: 
                print "Will consider the mean and width as parameter of interest"
                self.meanAsPOI = True
                self.widthAsPOI= True
            
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        if self.modelBuilder.out.var("r"):
            print "have r inside"
        else:
            self.modelBuilder.doVar("r[1,0,1000]")
        if self.muAsPOI: 
            poi = "r"
        if self.meanAsPOI: 
          if self.modelBuilder.out.var("mean_pole"):
              print "have mean inside"
          else:
              self.modelBuilder.doVar("mean_pole[500,0,600]")
              self.modelBuilder.doVar("sigma_pole[20,0,600]")
          poi = "mean_pole"
          poi += ",sigma_pole"
        
        self.modelBuilder.doSet("POI",poi)
        
HighmassModel= HighmassModel()
