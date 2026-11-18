import ROOT
import os
from math import sqrt
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class lepRecoSFProducerForFake(Module):
    def __init__(self):
        el_below_f = ["egammaEffi_ptBelow20.txt_EGM2D_UL2018.root"]
        el_above_f = ["egammaEffi_ptAbove20.txt_EGM2D_UL2018.root"]
        el_below_h = ["EGamma_SF2D"]
        el_above_h = ["EGamma_SF2D"]
        el_below_f = ["%s/src/TZWi/TopAnalysis/data/leptonRecoSF/" % os.environ['CMSSW_BASE'] + f for f in el_below_f]
        el_above_f = ["%s/src/TZWi/TopAnalysis/data/leptonRecoSF/" % os.environ['CMSSW_BASE'] + f for f in el_above_f]

        self.el_below_f = ROOT.std.vector(str)(len(el_below_f))
        self.el_above_f = ROOT.std.vector(str)(len(el_above_f))
        self.el_below_h = ROOT.std.vector(str)(len(el_below_f))
        self.el_above_h = ROOT.std.vector(str)(len(el_above_f))
        for i in range(len(el_below_f)):
            self.el_below_f[i] = el_below_f[i]
            self.el_below_h[i] = el_below_h[i]
        for i in range(len(el_above_f)):
            self.el_above_f[i] = el_above_f[i]
            self.el_above_h[i] = el_above_h[i]

        if "libTZWiTopAnalysis.so" not in ROOT.gSystem.GetLibraries():
            ROOT.gSystem.Load("libTZWiTopAnalysis.so")
        if not hasattr(ROOT, "LeptonEfficiencyCorrector"):
            ROOT.gROOT.ProcessLine(".L %s/src/TZWi/TopAnalysis/interface/LeptonEfficiencyCorrector.h" % os.environ['CMSSW_BASE'])

    def beginJob(self):
        self._worker_el_below = ROOT.LeptonEfficiencyCorrector(self.el_below_f,self.el_below_h)
        self._worker_el_above = ROOT.LeptonEfficiencyCorrector(self.el_above_f,self.el_above_h)

    def endJob(self):
        # Explicitly delete C++ worker objects to prevent ROOT GC segfaults at exit
        if hasattr(self, '_worker_el_below'): del self._worker_el_below
        if hasattr(self, '_worker_el_above'): del self._worker_el_above

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("FakeLepton_RecoSF", "F")
        self.out.branch("FakeLepton_RecoSFerr", "F")
        self.out.branch("FakeElectron_RecoSF", "F")
        self.out.branch("FakeElectron_RecoSFerr", "F")
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        sf_lepton1 = sf_lepton2 = sf_lepton3 = 1.
        sf_lepton1_err = sf_lepton2_err = sf_lepton3_err = 0.
        sf_el = 1.
        sf_el_err2 = 0.

        if abs(event._tree.b_out_Lepton1_pdgId) == 11:
            if event._tree.b_out_Lepton1_pt < 20:
                sf_lepton1 = self._worker_el_below.getSF(event._tree.b_out_Lepton1_pdgId, event._tree.b_out_Lepton1_pt, event._tree.b_out_Lepton1_eta)
                sf_lepton1_err = self._worker_el_below.getSFErr(event._tree.b_out_Lepton1_pdgId, event._tree.b_out_Lepton1_pt, event._tree.b_out_Lepton1_eta)
                sf_el *= sf_lepton1
                sf_el_err2 += sf_lepton1_err**2
            else:
                sf_lepton1 = self._worker_el_above.getSF(event._tree.b_out_Lepton1_pdgId, event._tree.b_out_Lepton1_pt, event._tree.b_out_Lepton1_eta)
                sf_lepton1_err = self._worker_el_above.getSFErr(event._tree.b_out_Lepton1_pdgId, event._tree.b_out_Lepton1_pt, event._tree.b_out_Lepton1_eta)
                sf_el *= sf_lepton1
                sf_el_err2 += sf_lepton1_err**2

        if abs(event._tree.b_out_Lepton2_pdgId) == 11:
            if event._tree.b_out_Lepton2_pt < 20:
                sf_lepton2 = self._worker_el_below.getSF(event._tree.b_out_Lepton2_pdgId, event._tree.b_out_Lepton2_pt, event._tree.b_out_Lepton2_eta)
                sf_lepton2_err = self._worker_el_below.getSFErr(event._tree.b_out_Lepton2_pdgId, event._tree.b_out_Lepton2_pt, event._tree.b_out_Lepton2_eta)
                sf_el *= sf_lepton2
                sf_el_err2 += sf_lepton2_err**2
            else:
                sf_lepton2 = self._worker_el_above.getSF(event._tree.b_out_Lepton2_pdgId, event._tree.b_out_Lepton2_pt, event._tree.b_out_Lepton2_eta)
                sf_lepton2_err = self._worker_el_above.getSFErr(event._tree.b_out_Lepton2_pdgId, event._tree.b_out_Lepton2_pt, event._tree.b_out_Lepton2_eta)
                sf_el *= sf_lepton2
                sf_el_err2 += sf_lepton2_err**2

        if abs(event._tree.b_out_Lepton3_pdgId) == 11:
            if event._tree.b_out_Lepton3_pt < 20:
                sf_lepton3 = self._worker_el_below.getSF(event._tree.b_out_Lepton3_pdgId, event._tree.b_out_Lepton3_pt, event._tree.b_out_Lepton3_eta)
                sf_lepton3_err = self._worker_el_below.getSFErr(event._tree.b_out_Lepton3_pdgId, event._tree.b_out_Lepton3_pt, event._tree.b_out_Lepton3_eta)
                sf_el *= sf_lepton3
                sf_el_err2 += sf_lepton3_err**2
            else:
                sf_lepton3 = self._worker_el_above.getSF(event._tree.b_out_Lepton3_pdgId, event._tree.b_out_Lepton3_pt, event._tree.b_out_Lepton3_eta)
                sf_lepton3_err = self._worker_el_above.getSFErr(event._tree.b_out_Lepton3_pdgId, event._tree.b_out_Lepton3_pt, event._tree.b_out_Lepton3_eta)
                sf_el *= sf_lepton3
                sf_el_err2 += sf_lepton3_err**2

        sf_lepton_total = sf_lepton1 * sf_lepton2 * sf_lepton3
        sf_lepton_err_total = sqrt( sf_lepton1_err**2 + sf_lepton2_err**2 + sf_lepton3_err**2 )
        sf_el_err = sqrt(sf_el_err2)

        self.out.fillBranch("FakeLepton_RecoSF", sf_lepton_total)
        self.out.fillBranch("FakeLepton_RecoSFerr", sf_lepton_err_total)
        self.out.fillBranch("FakeElectron_RecoSF", sf_el)
        self.out.fillBranch("FakeElectron_RecoSFerr", sf_el_err)
        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

lepRecoSF_forfake = lambda : lepRecoSFProducerForFake()
