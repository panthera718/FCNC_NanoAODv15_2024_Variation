import ROOT
import os
import math
import copy
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class JetVetoProducer(Module, object):

    def __init__(self, *args, **kwargs):
        self.jetVar = kwargs.get("jetVar")

        jveto_fN = ["hotjets-UL16.root"]
        jveto_h1 = ["h2hot_ul16_plus_hbm2_hbp12_qie11"]
        jveto_h2 = ["h2hot_mc"]

        jveto_f = ["%s/src/TZWi/TopAnalysis/data/JetVeto/" % os.environ['CMSSW_BASE'] + f for f in jveto_fN]

        self.jveto_f = ROOT.std.vector(str)(len(jveto_f))
        self.jveto_h1 = ROOT.std.vector(str)(len(jveto_f))
        self.jveto_h2 = ROOT.std.vector(str)(len(jveto_f))
        for i in range(len(jveto_f)): self.jveto_f[i] = jveto_f[i]; self.jveto_h1[i] = jveto_h1[i]; self.jveto_h2[i] = jveto_h2[i];

        if "libTZWiTopAnalysis.so" not in ROOT.gSystem.GetLibraries():
            ROOT.gSystem.Load("libTZWiTopAnalysis.so")
        if not hasattr(ROOT, "JetVetoCorrector"):
            ROOT.gROOT.ProcessLine(".L %s/src/TZWi/TopAnalysis/interface/JetVetoCorrector.h" % os.environ['CMSSW_BASE'])

    def beginJob(self):
        self._worker_jveto1 = ROOT.JetVetoCorrector(self.jveto_f,self.jveto_h1)
        self._worker_jveto2 = ROOT.JetVetoCorrector(self.jveto_f,self.jveto_h2)

    def endJob(self):
        # Explicitly delete C++ worker objects to prevent ROOT GC segfaults at exit
        if hasattr(self, '_worker_jveto1'): del self._worker_jveto1
        if hasattr(self, '_worker_jveto2'): del self._worker_jveto2

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        jetVar = self.jetVar

        self.out = wrappedOutputTree
#        self.out.branch("JetVeto", "F", lenVar="nGoodJet")
        self.out.branch("JetVetoSum_%sUp" % jetVar, "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        jetVar = self.jetVar

        jetveto = []
        jetvetosum = 0.

        for jindex in range(getattr(event._tree, "b_out_nGoodJet_%sUp" % jetVar)):
            jeta = (getattr(event._tree, "b_out_GoodJet_eta_%sUp" % jetVar))[jindex]
            jphi = (getattr(event._tree, "b_out_GoodJet_phi_%sUp" % jetVar))[jindex]
            jetvetoElement1 = self._worker_jveto1.getJV(jeta,jphi)
            jetvetoElement2 = self._worker_jveto2.getJV(jeta,jphi)
            jetveto.append(jetvetoElement1)
            jetveto.append(jetvetoElement2)
            jetvetosum += jetvetoElement1
            jetvetosum += jetvetoElement2

#        setattr(event._tree, "b_out_JetVeto", jetveto)
        setattr(event._tree, "b_out_JetVetoSum_%sUp" % jetVar, jetvetosum)

#        self.out.fillBranch("JetVeto", jetveto)
        self.out.fillBranch("JetVetoSum_%sUp" % jetVar, jetvetosum)

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

JVProducer_AbsoluteStat = lambda : JetVetoProducer(jetVar="AbsoluteStat")
JVProducer_AbsoluteScale = lambda : JetVetoProducer(jetVar="AbsoluteScale")
JVProducer_AbsoluteMPFBias = lambda : JetVetoProducer(jetVar="AbsoluteMPFBias")
JVProducer_Fragmentation = lambda : JetVetoProducer(jetVar="Fragmentation")
JVProducer_SinglePionECAL = lambda : JetVetoProducer(jetVar="SinglePionECAL")
JVProducer_SinglePionHCAL = lambda : JetVetoProducer(jetVar="SinglePionHCAL")
JVProducer_FlavorQCD = lambda : JetVetoProducer(jetVar="FlavorQCD")
JVProducer_FlavorPureQuark = lambda : JetVetoProducer(jetVar="FlavorPureQuark")
JVProducer_FlavorPureCharm = lambda : JetVetoProducer(jetVar="FlavorPureCharm")
JVProducer_FlavorPureBottom = lambda : JetVetoProducer(jetVar="FlavorPureBottom")
JVProducer_TimePtEta = lambda : JetVetoProducer(jetVar="TimePtEta")
JVProducer_RelativeJEREC1 = lambda : JetVetoProducer(jetVar="RelativeJEREC1")
JVProducer_RelativeJEREC2 = lambda : JetVetoProducer(jetVar="RelativeJEREC2")
JVProducer_RelativeJERHF = lambda : JetVetoProducer(jetVar="RelativeJERHF")
JVProducer_RelativePtBB = lambda : JetVetoProducer(jetVar="RelativePtBB")
JVProducer_RelativePtEC1 = lambda : JetVetoProducer(jetVar="RelativePtEC1")
JVProducer_RelativePtEC2 = lambda : JetVetoProducer(jetVar="RelativePtEC2")
JVProducer_RelativePtHF = lambda : JetVetoProducer(jetVar="RelativePtHF")
JVProducer_RelativeBal = lambda : JetVetoProducer(jetVar="RelativeBal")
JVProducer_RelativeSample = lambda : JetVetoProducer(jetVar="RelativeSample")
JVProducer_RelativeFSR = lambda : JetVetoProducer(jetVar="RelativeFSR")
JVProducer_RelativeStatFSR = lambda : JetVetoProducer(jetVar="RelativeStatFSR")
JVProducer_RelativeStatEC = lambda : JetVetoProducer(jetVar="RelativeStatEC")
JVProducer_RelativeStatHF = lambda : JetVetoProducer(jetVar="RelativeStatHF")
JVProducer_PileUpDataMC = lambda : JetVetoProducer(jetVar="PileUpDataMC")
JVProducer_PileUpPtRef = lambda : JetVetoProducer(jetVar="PileUpPtRef")
JVProducer_PileUpPtBB = lambda : JetVetoProducer(jetVar="PileUpPtBB")
JVProducer_PileUpPtEC1 = lambda : JetVetoProducer(jetVar="PileUpPtEC1")
JVProducer_PileUpPtEC2 = lambda : JetVetoProducer(jetVar="PileUpPtEC2")
JVProducer_PileUpPtHF = lambda : JetVetoProducer(jetVar="PileUpPtHF")
