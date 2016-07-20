import math
from HiggsAnalysis.CombinedLimit.PhysicsModel import *

### This is the base python class to study the SpinZero structure

class SpinZeroHiggs(PhysicsModel):
    def __init__(self):
        self.mHRange = []
        self.muAsPOI = False
        self.muFloating = False
        self.fg4fixed = False
        self.fg2POI = False
        self.fg4POI = True 
        self.fcpPOI = False
        self.phiPOI = False
        self.cPOI = False
        self.phiFloating = False
        self.phi2Floating = False
        self.poiMap = []
        self.pois = {}
        self.verbose = False
    def setModelBuilder(self, modelBuilder):
        PhysicsModel.setModelBuilder(self,modelBuilder)
        self.modelBuilder.doModelBOnly = False

    def getYieldScale(self,bin,process):
        "Split in production and decay, and call getHiggsSignalYieldScale; return 1 for backgrounds "
        #print "Bin ",bin
        #print "Process ",process
        if self.DC.isSignal[process]:
            if self.muFloating:
                self.my_norm = "r"
            else:
                self.my_norm = 1 
        
            print "Process {0} will scale by {1}".format(process,self.my_norm)
            return self.my_norm
        
        elif not self.DC.isSignal[process]: return 1
            

    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            #if po == "muAsPOI":
            if 'muAsPOI' in po: 
                print "Will consider the signal strength as a parameter of interest"
                self.muAsPOI = True
                self.muFloating = True
            #if po == "muFloating":
            if 'muFloating' in po: 
                print "Will consider the signal strength as a floating parameter (as a paraeter of interest if --PO muAsPOI is specified, as a nuisance otherwise)"
                self.muFloating = True
            #if po == "fg4fixed":
            if 'fg4fixed' in po: 
                print "Will fix CMS_zz4l_fg4 to 0 and make signal strength as a parameter of interest"
                self.muAsPOI = True
                self.muFloating = True
                self.fg4fixed = True
            if '2D' in po: 
                print "Will set CMS_zz4l_fg2 and CMS_zz4l_fg4 as parameter of interest"
                self.fg4fixed = False 
                self.fg2POI = True 
                self.phi2Floating = True 
                self.phiFloating = True 
            if 'fcp' in po: 
                print "Will set fcp as parameter of interest"
                self.fg4fixed = False 
                self.fg4POI = False 
                self.fcpPOI = True 
            #if po == "muFloating,phiAsPOI":
            if 'phiAsPOI' in po: 
                print "Will consider the phase as a parameter of interest"
                self.phiPOI = True
                self.phiFloating = True
            if 'phiFloating' in po: 
                print "Will consider the phase as a floating parameter"
                self.phiFloating = True
            if 'cPOI' in po: 
                print "Will consider the czz_ww as a parameter of interest"
                self.cPOI = True
            
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        if self.fg4fixed:
            self.modelBuilder.doVar("CMS_zz4l_fg4[0]")
            self.modelBuilder.doVar("r[1,0,4]")
            print "Fixing CMS_zz4l_fg4"
            poi = "r"
        else:
            if self.modelBuilder.out.var("CMS_zz4l_fg4"):
                print "have fg4 inside"
            else:
                self.modelBuilder.doVar("CMS_zz4l_fg4[0.,-1,1]")
            if self.fg4POI:
                poi = "CMS_zz4l_fg4"
            if self.cPOI:
                if self.modelBuilder.out.var("cww_zz"):
                   print "have czz_ww inside"
                else:
                   self.modelBuilder.doVar("cww_zz[0.5,-10,10]")
                poi += ",cww_zz"

            if self.fg2POI:
               if self.modelBuilder.out.var("CMS_zz4l_fg2"):
                   print "have fg2 inside"
               else:
                   self.modelBuilder.doVar("CMS_zz4l_fg2[0.,0,1]")
               poi += ",CMS_zz4l_fg2"
            if self.fcpPOI:
               if self.modelBuilder.out.var("fcp"):
                   print "have fcp inside"
               else:
                   self.modelBuilder.doVar("fcp[0.,-1,1]")
               poi = "fcp"
            if self.muFloating:
                self.modelBuilder.doVar("r[1,0,2000]")
                if self.muAsPOI:
                    print "Treating r as a POI"
                    poi += ",r"
                else:
                    self.modelBuilder.out.var("r").setAttribute("flatParam")
            if self.phiFloating:
                if self.modelBuilder.out.var("CMS_zz4l_fg4phi"):
                    print "have fg4phi inside"
                else: 
                    self.modelBuilder.doVar("CMS_zz4l_fg4phi[0.,-3.1415926,3.1415926]")
                if self.phiPOI:
                    poi += ",CMS_zz4l_fg4phi"
                else:
                    self.modelBuilder.out.var("CMS_zz4l_fg4phi").setAttribute("flatParam")
            if self.phi2Floating:
                #self.modelBuilder.doVar("CMS_zz4l_fg4phi[0.,-math.pi,math.pi]")
                if self.modelBuilder.out.var("CMS_zz4l_fg2phi"):
                    print "have fg2phi inside"
                else: 
                    self.modelBuilder.doVar("CMS_zz4l_fg2phi[0.,-3.1415926,3.1415926]")
                self.modelBuilder.out.var("CMS_zz4l_fg2phi").setAttribute("flatParam")
        
        self.modelBuilder.doSet("POI",poi)
        
spinZeroHiggs = SpinZeroHiggs()
