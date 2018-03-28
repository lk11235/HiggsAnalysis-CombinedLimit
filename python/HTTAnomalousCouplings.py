#conventions:
#sigma1 is for a1=1
#sigma3 is for a3=1 (similarly for a2, L1, L1Zg)
#sigmaa1a3int is for the mixture sample when a1=1, a3=(our convention for fa3=0.5 samples)

xsecs = {
    "sigma1_HZZ": 290.58626,
    "sigma3_HZZ": 44.670158,
    "sigma1_VBF": 968.674,
    "sigma3_VBF": 10909.54,
    "sigma1_ZH":  9022.36,
    "sigma3_ZH":  434763.7,
    "sigma1_WH":  30998.54,
    "sigma3_WH":  2028656,
    "sigma2_HZZ": 105.85594,
    "sigmaL1_HZZ": 1.9846071e-06,
    "sigma2_VBF": 13102.71,
    "sigmaL1_VBF": 2.08309E-4,
    "sigma2_ZH": 713123,
    "sigmaL1_ZH": 33652.46e-6,
    "sigma2_WH": 3106339,
    "sigmaL1_WH": 11234.91e-5,
    "sigmaa1a3int_VBF": 1937.15,
    "sigmaa1a3int_ZH": 18044.72,
    "sigmaa1a3int_WH": 61997.07,
    # for ggH what is listed are yields in histos (all normliazed to powheg xsection):
    "sigmaa1a2int_VBF": 2207.73 ,
    "sigmaa1a2int_ZH": 4258.966,
    "sigmaa1a2int_WH": 16486.68,
    "sigmaa1L1int_VBF": 2861.21349769,
    "sigmaa1L1int_ZH": 6852.307, 
    "sigmaa1L1int_WH": 25302.37, 
    "sigmaa1L1Zgint_VBF": 1410.5494,
    "sigmaa1L1Zgint_ZH": 21952.490480117736,
    "sigmaa1L1Zgint_WH": 30998.54,     #this is the same as SM, because L1Zg doesn't contribute to WH
    "yield_Powheg_ggH": 5.127128e+00,
    "yield_SM_ggH": 7.937751e+01,
    "yield_BSM_ggH": 8.952848e+01,
    "sigma_Powheg_ggH": 48.58,
    "sigma_SM_ggH": 15980,
    "sigma_BSM_ggH": 15981,
    "BR_H_tt": 0.0627
}

class Anomalous_Interference_JHU_rw(PhysicsModelBase_NiceSubclasses):
    def __init__(self, *args, **kwargs):
        self.anomalouscoupling = None
        self.dofa3gg = None
        super(Anomalous_Interference_JHU_rw, self).__init__(*args, **kwargs)

    def getPOIList(self):
        myxsecs = {
            "sigma1_HZZ": xsecs["sigma1_HZZ"]
            "sigma1_VBF": xsecs["sigma1_VBF"]
            "sigma1_ZH": xsecs["sigma1_ZH"]
            "sigma1_WH": xsecs["sigma1_WH"]
        }
        if self.anomalouscoupling == "fa3":
            myxsecs.update({
                "sigmai_HZZ": xsecs["sigma3_HZZ"]
                "sigmai_VBF": xsecs["sigma3_VBF"]
                "sigmai_ZH": xsecs["sigma3_ZH"]
                "sigmai_WH": xsecs["sigma3_WH"]
                "sigmaint_HZZ": xsecs["sigmaa1a3int_HZZ"]
                "sigmaint_VBF": xsecs["sigmaa1a3int_VBF"]
                "sigmaint_ZH": xsecs["sigmaa1a3int_ZH"]
                "sigmaint_WH": xsecs["sigmaa1a3int_WH"]
            })
        elif self.anomalouscoupling == "fa2":
            myxsecs.update({
                "sigmai_HZZ": xsecs["sigma2_HZZ"]
                "sigmai_VBF": xsecs["sigma2_VBF"]
                "sigmai_ZH": xsecs["sigma2_ZH"]
                "sigmai_WH": xsecs["sigma2_WH"]
                "sigmaint_HZZ": xsecs["sigmaa1a2int_HZZ"]
                "sigmaint_VBF": xsecs["sigmaa1a2int_VBF"]
                "sigmaint_ZH": xsecs["sigmaa1a2int_ZH"]
                "sigmaint_WH": xsecs["sigmaa1a2int_WH"]
            })
        elif self.anomalouscoupling == "fL1":
            myxsecs.update({
                "sigmai_HZZ": xsecs["sigmaL1_HZZ"]
                "sigmai_VBF": xsecs["sigmaL1_VBF"]
                "sigmai_ZH": xsecs["sigmaL1_ZH"]
                "sigmai_WH": xsecs["sigmaL1_WH"]
                "sigmaint_HZZ": xsecs["sigmaa1L1int_HZZ"]
                "sigmaint_VBF": xsecs["sigmaa1L1int_VBF"]
                "sigmaint_ZH": xsecs["sigmaa1L1int_ZH"]
                "sigmaint_WH": xsecs["sigmaa1L1int_WH"]
            })
        elif self.anomalouscoupling == "fL1Zg":
            myxsecs.update({
                "sigmai_HZZ": xsecs["sigmaL1Zg_HZZ"]
                "sigmai_VBF": xsecs["sigmaL1Zg_VBF"]
                "sigmai_ZH": xsecs["sigmaL1Zg_ZH"]
                "sigmai_WH": xsecs["sigmaL1Zg_WH"]
                "sigmaint_HZZ": xsecs["sigmaa1L1Zgint_HZZ"]
                "sigmaint_VBF": xsecs["sigmaa1L1Zgint_VBF"]
                "sigmaint_ZH": xsecs["sigmaa1L1Zgint_ZH"]
                "sigmaint_WH": xsecs["sigmaa1L1Zgint_WH"]
            })

        """Create POI and other parameters, and define the POI set."""
        if not self.modelBuilder.out.var("CMS_zz4l_fai1"):
            self.modelBuilder.doVar("CMS_zz4l_fai1[0.0,-1.0,1.0]")
        if not self.modelBuilder.out.var("RF"):
            self.modelBuilder.doVar("RF[1.0,0,10]")
        if not self.modelBuilder.out.var("muTT"):
            self.modelBuilder.doVar("muTT[1.0,0,10]")
        if not self.modelBuilder.out.var("RV"):
            self.modelBuilder.doVar("RV[1.0,0.0,10.0]")

        self.modelBuilder.doVar('expr::a1("sqrt(1-abs(@0))", CMS_zz4l_fai1)')
        self.modelBuilder.doVar('expr::ai("(@0>0 ? 1 : -1) * sqrt(abs(@0)*{sigma1_HZZ}/{sigmai_HZZ})", CMS_zz4l_fai1)'.format(**myxsecs))

        self.modelBuilder.factory_('expr::smCoupling_VBF("@0*@3*@1**2 - @0*@1*@2*@3*sqrt({sigmai_VBF}/{sigma1_VBF})", RV,a1,ai,muTT)'.format(**xsecs))
        self.modelBuilder.factory_('expr::smCoupling_ZH("@0*@3*@1**2 - @0*@1*@2*@3*sqrt({sigmai_ZH}/{sigma1_ZH})", RV,a1,ai,muTT)'.format(**xsecs))
        self.modelBuilder.factory_('expr::smCoupling_WH("@0*@3*@1**2 - @0*@1*@2*@3*sqrt({sigmai_WH}/{sigma1_WH})", RV,a1,ai,muTT)'.format(**xsecs))

        self.modelBuilder.factory_('expr::bsmCoupling_VBF("@0*@3*@1**2*{sigmai_VBF}/{sigma1_VBF} - @0*@1*@2*@3*sqrt({sigmai_VBF}/{sigma1_VBF})", RV,ai,a1,muTT)'.format(**xsecs))
        self.modelBuilder.factory_('expr::bsmCoupling_ZH("@0*@3*@1**2*{sigmai_ZH}/{sigma1_ZH} - @0*@1*@2*@3*sqrt({sigmai_ZH}/{sigma1_ZH})", RV,ai,a1,muTT)'.format(**xsecs))
        self.modelBuilder.factory_('expr::bsmCoupling_WH("@0*@3*@1**2*{sigmai_WH}/{sigma1_WH} - @0*@1*@2*@3*sqrt({sigmai_WH}/{sigma1_WH})", RV,ai,a1,muTT)'.format(**xsecs))

        self.modelBuilder.factory_('expr::intCoupling_VBF("@0*@1*@2*@3*sqrt({sigmai_VBF}/{sigma1_VBF})*{sigmaint_VBF}/{sigma1_VBF}", RV,a1,ai,muTT)'.format(**xsecs))
        self.modelBuilder.factory_('expr::intCoupling_ZH("@0*@1*@2*@3*sqrt({sigmai_ZH}/{sigma1_ZH})*{sigmaint_ZH}/{sigma1_ZH}", RV,a1,ai,muTT)'.format(**xsecs))
        self.modelBuilder.factory_('expr::intCoupling_WH("@0*@1*@2*@3*sqrt({sigmai_WH}/{sigma1_WH})*{sigmaint_WH}/{sigma1_WH}", RV,a1,ai,muTT)'.format(**xsecs))

        pois = ["CMS_zz4l_fai1","RV","RF","muTT"]

        if self.dofa3gg:
            if not self.modelBuilder.out.var("fa3_ggH"):
                self.modelBuilder.doVar("fa3_ggH[0.0,0.0,1.0]")
            self.modelBuilder.doVar('expr::a1_ggH("sqrt(1-abs(@0))", fa3_ggH)')
            self.modelBuilder.doVar('expr::ai_ggH("sqrt(abs(@0))", fa3_ggH)'.format(**xsecs))
            self.modelBuilder.factory_('expr::smCoupling_ggH("@0*@2*@1**2", RF,a1_ggH,muTT)'.format(**xsecs))
            self.modelBuilder.factory_('expr::bsmCoupling_ggH("@0*@2*@1**2", RF,a3_ggH,muTT)'.format(**xsecs))
            pois.append("fa3_ggH")
        else:
            self.modelBuilder.factory_('expr::smCoupling_ggH("@0*@2", RF,muTT)'.format(**xsecs))
            

        return pois + super(Anomalous_Interference_JHU_rw, self).getPOIList()

    def getYieldScale(self,bin,process):
        if process in ["GGH2Jets_sm_M",]:
            return 'smCoupling_ggH'
        if process in ["GGH2Jets_pseudoscalar_M",]:
            if not self.dofa3gg: raise ValueError("You have GGH2Jets_pseudoscalar_M, but did not set dofa3gg=true")
            return 'bsmCoupling_ggH'
        if process in ["reweighted_qqH_htt_0PM",]:
            return 'smCoupling_VBF'
        if process in ["reweighted_WH_htt_0PM",]:
            return 'smCoupling_WH'
        if process in ["reweighted_ZH_htt_0PM",]:
            return 'smCoupling_ZH'
        if process in ["reweighted_qqH_htt_0L1",]:
            return 'bsmCoupling_VBF'
        if process in ["reweighted_WH_htt_0L1",]:
            return 'bsmCoupling_WH'
        if process in ["reweighted_ZH_htt_0L1",]:
            return 'bsmCoupling_ZH'
        if process in ["reweighted_qqH_htt_0L1f05ph0"]:
            return 'intCoupling_VBF'
        if process in ["reweighted_ZH_htt_0L1f05ph0"]:
            return 'intCoupling_ZH'
        if process in ["reweighted_WH_htt_0L1f05ph0"]:
            return 'intCoupling_WH'
        return super(Anomalous_Interference_JHU_rw, self).getYieldScale(bin, process)

    def processPhysicsOptions(self,physOptions):
        processed = super(Anomalous_Interference_JHU_rw, self).processPhysicsOptions(physOptions)
        for po in physOptions:
            if po in ("fa3", "fa2", "fL1", "fL1Zg"):
                if self.anomalouscoupling is not None:
                    raise ValueError("Multiple anomalous couplings provided as physics options ({}, {})".format(self.anomalouscoupling, po))
                self.anomalouscoupling = po
                processed.append(po)
            if po.lower() == "dofa3gg=true":
                self.dofa3gg = True
                processed.append(po)
            if po.lower() == "dofa3gg=false":
                self.dofa3gg = False
                processed.append(po)

        if self.anomalouscoupling is None:
            raise ValueError("Have to provide an anomalous coupling as a physics option.  Choices: fa3, fa2, fL1, fL1Zg")
        if self.dofa3gg is None:  #not set
            dofa3gg = (self.anomalouscoupling == "fa3")

        return processed

class Anomalous_Interference_JHU_rw_HTTHZZ(Anomalous_Interference_JHU_rw, MultiSignalSpinZeroHiggs): pass

Anomalous_Interference_JHU_rw = Anomalous_Interference_JHU_rw()
Anomalous_Interference_JHU_rw_HTTHZZ = Anomalous_Interference_JHU_rw_HTTHZZ()
