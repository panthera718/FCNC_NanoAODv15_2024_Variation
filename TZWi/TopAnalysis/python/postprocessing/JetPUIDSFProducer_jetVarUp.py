import ROOT
import os
import math
import copy
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class JetPUIDSFProducer(Module, object):

    def __init__(self, *args, **kwargs):

        self.jetVar = kwargs.get("jetVar")

        jpuidsf_fN = ["PUID_106XTraining_ULRun2_EffSFandUncties_v1.root"]
        jpuidsf_h = ["h2_eff_sfUL2018_L"]
        jpuidsf_herr = ["h2_eff_sfUL2018_L_Systuncty"]

        jpuidsf_f = ["%s/src/TZWi/TopAnalysis/data/JetPUIDSF/" % os.environ['CMSSW_BASE'] + f for f in jpuidsf_fN]

        self.jpuidsf_f = ROOT.std.vector(str)(len(jpuidsf_f))
        self.jpuidsf_h = ROOT.std.vector(str)(len(jpuidsf_f))
        self.jpuidsf_herr = ROOT.std.vector(str)(len(jpuidsf_f))
        for i in range(len(jpuidsf_f)): self.jpuidsf_f[i] = jpuidsf_f[i]; self.jpuidsf_h[i] = jpuidsf_h[i]; self.jpuidsf_herr[i] = jpuidsf_herr[i];

        if "/JetPUIDSFCorrector_cc.so" not in ROOT.gSystem.GetLibraries():
            print("Load C++ Worker")
            ROOT.gROOT.ProcessLine(".L %s/src/TZWi/TopAnalysis/python/postprocessing/helpers/JetPUIDSFCorrector.cc+" % os.environ['CMSSW_BASE'])

    def beginJob(self):
        self._worker_jpuidsf = ROOT.JetPUIDSFCorrector(self.jpuidsf_f,self.jpuidsf_h)
        self._worker_jpuidsferr = ROOT.JetPUIDSFCorrector(self.jpuidsf_f,self.jpuidsf_herr)

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        jetVar = self.jetVar
#        self.out.branch("JetPUID_eachSF", "F", lenVar="nGoodJet")
#        self.out.branch("JetPUID_eachSFerr", "F", lenVar="nGoodJet")
        self.out.branch("JetPUID_SF_%sUp" % jetVar, "F")
        self.out.branch("JetPUID_SFerr_%sUp" % jetVar, "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jetVar = self.jetVar

        jetpuid_eachsf = []
        jetpuid_eachsferr = []
        jetpuid_sf = 1.
        jetpuid_sferr_sqsum = 0.

        for jetindex in range(event._tree.b_out_nGoodJet):
            jpt = (event._tree.b_out_GoodJet_pt)[jetindex]
            jeta = (event._tree.b_out_GoodJet_eta)[jetindex]
            if jpt < 50.:
                jetpuidsfElement = self._worker_jpuidsf.getJPUIDSF(jpt,jeta)
                jetpuidsferrElement = self._worker_jpuidsferr.getJPUIDSF(jpt,jeta)
                jetpuid_eachsf.append(jetpuidsfElement)
                jetpuid_eachsferr.append(jetpuidsferrElement)
                jetpuid_sf *= jetpuidsfElement
                jetpuid_sferr_sqsum += jetpuidsferrElement**2

        jetpuid_sferr = math.sqrt(jetpuid_sferr_sqsum)

#        setattr(event._tree, "b_out_JetPUID_eachSF", jetpuid_eachsf)
#        setattr(event._tree, "b_out_JetPUID_eachSFerr", jetpuid_eachsferr)
        setattr(event._tree, "b_out_JetPUID_SF_%sUp" % jetVar, jetpuid_sf)
        setattr(event._tree, "b_out_JetPUID_SFerr_%sUp" % jetVar, jetpuid_sferr)

#        self.out.fillBranch("JetPUID_eachSF", jetpuid_eachsf)
#        self.out.fillBranch("JetPUID_eachSFerr", jetpuid_eachsferr)
        self.out.fillBranch("JetPUID_SF_%sUp" % jetVar, jetpuid_sf)
        self.out.fillBranch("JetPUID_SFerr_%sUp" % jetVar, jetpuid_sferr)

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

JPUIDSFProducer_AbsoluteStat = lambda : JetPUIDSFProducer(jetVar="AbsoluteStat")
JPUIDSFProducer_AbsoluteScale = lambda : JetPUIDSFProducer(jetVar="AbsoluteScale")
JPUIDSFProducer_AbsoluteMPFBias = lambda : JetPUIDSFProducer(jetVar="AbsoluteMPFBias")
JPUIDSFProducer_Fragmentation = lambda : JetPUIDSFProducer(jetVar="Fragmentation")
JPUIDSFProducer_SinglePionECAL = lambda : JetPUIDSFProducer(jetVar="SinglePionECAL")
JPUIDSFProducer_SinglePionHCAL = lambda : JetPUIDSFProducer(jetVar="SinglePionHCAL")
JPUIDSFProducer_FlavorPureGluon = lambda : JetPUIDSFProducer(jetVar="FlavorPureGluon")
JPUIDSFProducer_FlavorPureQuark = lambda : JetPUIDSFProducer(jetVar="FlavorPureQuark")
JPUIDSFProducer_FlavorPureCharm = lambda : JetPUIDSFProducer(jetVar="FlavorPureCharm")
JPUIDSFProducer_FlavorPureBottom = lambda : JetPUIDSFProducer(jetVar="FlavorPureBottom")
JPUIDSFProducer_TimePtEta = lambda : JetPUIDSFProducer(jetVar="TimePtEta")
JPUIDSFProducer_RelativeJEREC1 = lambda : JetPUIDSFProducer(jetVar="RelativeJEREC1")
JPUIDSFProducer_RelativeJEREC2 = lambda : JetPUIDSFProducer(jetVar="RelativeJEREC2")
JPUIDSFProducer_RelativeJERHF = lambda : JetPUIDSFProducer(jetVar="RelativeJERHF")
JPUIDSFProducer_RelativePtBB = lambda : JetPUIDSFProducer(jetVar="RelativePtBB")
JPUIDSFProducer_RelativePtEC1 = lambda : JetPUIDSFProducer(jetVar="RelativePtEC1")
JPUIDSFProducer_RelativePtEC2 = lambda : JetPUIDSFProducer(jetVar="RelativePtEC2")
JPUIDSFProducer_RelativePtHF = lambda : JetPUIDSFProducer(jetVar="RelativePtHF")
JPUIDSFProducer_RelativeBal = lambda : JetPUIDSFProducer(jetVar="RelativeBal")
JPUIDSFProducer_RelativeSample = lambda : JetPUIDSFProducer(jetVar="RelativeSample")
JPUIDSFProducer_RelativeFSR = lambda : JetPUIDSFProducer(jetVar="RelativeFSR")
JPUIDSFProducer_RelativeStatFSR = lambda : JetPUIDSFProducer(jetVar="RelativeStatFSR")
JPUIDSFProducer_RelativeStatEC = lambda : JetPUIDSFProducer(jetVar="RelativeStatEC")
JPUIDSFProducer_RelativeStatHF = lambda : JetPUIDSFProducer(jetVar="RelativeStatHF")
JPUIDSFProducer_PileUpDataMC = lambda : JetPUIDSFProducer(jetVar="PileUpDataMC")
JPUIDSFProducer_PileUpPtRef = lambda : JetPUIDSFProducer(jetVar="PileUpPtRef")
JPUIDSFProducer_PileUpPtBB = lambda : JetPUIDSFProducer(jetVar="PileUpPtBB")
JPUIDSFProducer_PileUpPtEC1 = lambda : JetPUIDSFProducer(jetVar="PileUpPtEC1")
JPUIDSFProducer_PileUpPtEC2 = lambda : JetPUIDSFProducer(jetVar="PileUpPtEC2")
JPUIDSFProducer_PileUpPtHF = lambda : JetPUIDSFProducer(jetVar="PileUpPtHF")
