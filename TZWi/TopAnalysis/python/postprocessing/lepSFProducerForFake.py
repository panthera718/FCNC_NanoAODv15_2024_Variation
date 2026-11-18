import ROOT
import os
from math import sqrt
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class lepSFProducerForFake(Module):
    def __init__(self):
        mu_f = ["NUM_LeptonMvaMedium_DEN_TrackerMuons_abseta_pt.root"]
        mu_h = ["NUM_LeptonMvaMedium_DEN_TrackerMuons_abseta_pt"]
        mu_h_stat = ["NUM_LeptonMvaMedium_DEN_TrackerMuons_abseta_pt_stat"]
        mu_h_syst = ["NUM_LeptonMvaMedium_DEN_TrackerMuons_abseta_pt_syst"]
        el_f = ["egammaEffi.txt_EGM2D.root"]
        el_h = ["EGamma_SF2D"]
        el_h_statErr = ["stat"]
        el_h_systErr = ["sys"]
        mu_f = ["%s/src/TZWi/TopAnalysis/data/leptonSF/" % os.environ['CMSSW_BASE'] + f for f in mu_f]
        el_f = ["%s/src/TZWi/TopAnalysis/data/leptonSF/" % os.environ['CMSSW_BASE'] + f for f in el_f]

        self.el_f = ROOT.std.vector(str)(len(el_f))
        self.el_h = ROOT.std.vector(str)(len(el_f))
        self.el_h_statErr = ROOT.std.vector(str)(len(el_f))
        self.el_h_systErr = ROOT.std.vector(str)(len(el_f))
        for i in range(len(el_f)):
            self.el_f[i] = el_f[i]
            self.el_h[i] = el_h[i]
            self.el_h_statErr[i] = el_h_statErr[i]
            self.el_h_systErr[i] = el_h_systErr[i]
        self.mu_f = ROOT.std.vector(str)(len(mu_f))
        self.mu_h = ROOT.std.vector(str)(len(mu_f))
        self.mu_h_stat = ROOT.std.vector(str)(len(mu_f))
        self.mu_h_syst = ROOT.std.vector(str)(len(mu_f))
        for i in range(len(mu_f)):
            self.mu_f[i] = mu_f[i]
            self.mu_h[i] = mu_h[i]
            self.mu_h_stat[i] = mu_h_stat[i]
            self.mu_h_syst[i] = mu_h_syst[i]

        if "libTZWiTopAnalysis.so" not in ROOT.gSystem.GetLibraries():
            ROOT.gSystem.Load("libTZWiTopAnalysis.so")
        if not hasattr(ROOT, "LeptonEfficiencyCorrector"):
            ROOT.gROOT.ProcessLine(".L %s/src/TZWi/TopAnalysis/interface/LeptonEfficiencyCorrector.h" % os.environ['CMSSW_BASE'])

    def beginJob(self):
        self._worker_el = ROOT.LeptonEfficiencyCorrector(self.el_f,self.el_h)
        self._worker_el_StatErr = ROOT.LeptonEfficiencyCorrector(self.el_f,self.el_h_statErr)
        self._worker_el_SystErr = ROOT.LeptonEfficiencyCorrector(self.el_f,self.el_h_systErr)
        self._worker_mu = ROOT.LeptonEfficiencyCorrector(self.mu_f,self.mu_h)
        self._worker_mu_Stat = ROOT.LeptonEfficiencyCorrector(self.mu_f,self.mu_h_stat)
        self._worker_mu_Syst = ROOT.LeptonEfficiencyCorrector(self.mu_f,self.mu_h_syst)

    def endJob(self):
        # Explicitly delete C++ worker objects to prevent ROOT GC segfaults at exit
        if hasattr(self, '_worker_el'): del self._worker_el
        if hasattr(self, '_worker_el_StatErr'): del self._worker_el_StatErr
        if hasattr(self, '_worker_el_SystErr'): del self._worker_el_SystErr
        if hasattr(self, '_worker_mu'): del self._worker_mu
        if hasattr(self, '_worker_mu_Stat'): del self._worker_mu_Stat
        if hasattr(self, '_worker_mu_Syst'): del self._worker_mu_Syst

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("FakeLepton_SF", "F")
        self.out.branch("FakeLepton_SFerr", "F")
        self.out.branch("FakeLepton_SFerrStat", "F")
        self.out.branch("FakeLepton_SFerrSyst", "F")
        #self.out.branch("FakeLepton1_SF", "F")
        #self.out.branch("FakeLepton1_SFerr", "F")
        #self.out.branch("FakeLepton2_SF", "F")
        #self.out.branch("FakeLepton2_SFerr", "F")
        #self.out.branch("FakeLepton3_SF", "F")
        #self.out.branch("FakeLepton3_SFerr", "F")
        self.out.branch("FakeElectron_SF", "F")
        self.out.branch("FakeElectron_SFerr", "F")
        self.out.branch("FakeElectron_SFerrStat", "F")
        self.out.branch("FakeElectron_SFerrSyst", "F")
        self.out.branch("FakeMuon_SF", "F")
        self.out.branch("FakeMuon_SFerr", "F")
        self.out.branch("FakeMuon_SFerrStat", "F")
        self.out.branch("FakeMuon_SFerrSyst", "F")
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        sf_lepton1 = sf_lepton2 = sf_lepton3 = 1.
        sf_lepton1_err = sf_lepton2_err = sf_lepton3_err = 0.
        sf_lepton1_errStat = sf_lepton2_errStat = sf_lepton3_errStat = 0.
        sf_lepton1_errSyst = sf_lepton2_errSyst = sf_lepton3_errSyst = 0.
        sf_el = sf_mu = 1.
        sf_el_err2 = sf_mu_err2 = 0.
        sf_el_err2Stat = sf_mu_err2Stat = 0.
        sf_el_err2Syst = sf_mu_err2Syst = 0.

        if abs(event._tree.b_out_Lepton1_pdgId) == 11:
            sf_lepton1 = self._worker_el.getSF(event._tree.b_out_Lepton1_pdgId, event._tree.b_out_Lepton1_pt, event._tree.b_out_Lepton1_eta)
            sf_lepton1_err = self._worker_el.getSFErr(event._tree.b_out_Lepton1_pdgId, event._tree.b_out_Lepton1_pt, event._tree.b_out_Lepton1_eta)
            sf_lepton1_errStat = self._worker_el_StatErr.getSF(event._tree.b_out_Lepton1_pdgId, event._tree.b_out_Lepton1_pt, event._tree.b_out_Lepton1_eta)
            sf_lepton1_errSyst = self._worker_el_SystErr.getSF(event._tree.b_out_Lepton1_pdgId, event._tree.b_out_Lepton1_pt, event._tree.b_out_Lepton1_eta)
            sf_el *= sf_lepton1
            sf_el_err2 += sf_lepton1_err**2
            sf_el_err2Stat += sf_lepton1_errStat**2
            sf_el_err2Syst += sf_lepton1_errSyst**2
        elif abs(event._tree.b_out_Lepton1_pdgId) == 13:
            sf_lepton1 = self._worker_mu.getSF(event._tree.b_out_Lepton1_pdgId, event._tree.b_out_Lepton1_pt, event._tree.b_out_Lepton1_eta)
            sf_lepton1_err = self._worker_mu.getSFErr(event._tree.b_out_Lepton1_pdgId, event._tree.b_out_Lepton1_pt, event._tree.b_out_Lepton1_eta)
            sf_lepton1_errStat = self._worker_mu_Stat.getSFErr(event._tree.b_out_Lepton1_pdgId, event._tree.b_out_Lepton1_pt, event._tree.b_out_Lepton1_eta)
            sf_lepton1_errSyst = self._worker_mu_Syst.getSFErr(event._tree.b_out_Lepton1_pdgId, event._tree.b_out_Lepton1_pt, event._tree.b_out_Lepton1_eta)
            sf_mu *= sf_lepton1
            sf_mu_err2 += sf_lepton1_err**2
            sf_mu_err2Stat += sf_lepton1_errStat**2
            sf_mu_err2Syst += sf_lepton1_errSyst**2

        if abs(event._tree.b_out_Lepton2_pdgId) == 11:
            sf_lepton2 = self._worker_el.getSF(event._tree.b_out_Lepton2_pdgId, event._tree.b_out_Lepton2_pt, event._tree.b_out_Lepton2_eta)
            sf_lepton2_err = self._worker_el.getSFErr(event._tree.b_out_Lepton2_pdgId, event._tree.b_out_Lepton2_pt, event._tree.b_out_Lepton2_eta)
            sf_lepton2_errStat = self._worker_el_StatErr.getSF(event._tree.b_out_Lepton2_pdgId, event._tree.b_out_Lepton2_pt, event._tree.b_out_Lepton2_eta)
            sf_lepton2_errSyst = self._worker_el_SystErr.getSF(event._tree.b_out_Lepton2_pdgId, event._tree.b_out_Lepton2_pt, event._tree.b_out_Lepton2_eta)
            sf_el *= sf_lepton2
            sf_el_err2 += sf_lepton2_err**2
            sf_el_err2Stat += sf_lepton2_errStat**2
            sf_el_err2Syst += sf_lepton2_errSyst**2
        elif abs(event._tree.b_out_Lepton2_pdgId) == 13:
            sf_lepton2 = self._worker_mu.getSF(event._tree.b_out_Lepton2_pdgId, event._tree.b_out_Lepton2_pt, event._tree.b_out_Lepton2_eta)
            sf_lepton2_err = self._worker_mu.getSFErr(event._tree.b_out_Lepton2_pdgId, event._tree.b_out_Lepton2_pt, event._tree.b_out_Lepton2_eta)
            sf_lepton2_errStat = self._worker_mu_Stat.getSFErr(event._tree.b_out_Lepton2_pdgId, event._tree.b_out_Lepton2_pt, event._tree.b_out_Lepton2_eta)
            sf_lepton2_errSyst = self._worker_mu_Syst.getSFErr(event._tree.b_out_Lepton2_pdgId, event._tree.b_out_Lepton2_pt, event._tree.b_out_Lepton2_eta)
            sf_mu *= sf_lepton2
            sf_mu_err2 += sf_lepton2_err**2
            sf_mu_err2Stat += sf_lepton2_errStat**2
            sf_mu_err2Syst += sf_lepton2_errSyst**2

        if abs(event._tree.b_out_Lepton3_pdgId) == 11:
            sf_lepton3 = self._worker_el.getSF(event._tree.b_out_Lepton3_pdgId, event._tree.b_out_Lepton3_pt, event._tree.b_out_Lepton3_eta)
            sf_lepton3_err = self._worker_el.getSFErr(event._tree.b_out_Lepton3_pdgId, event._tree.b_out_Lepton3_pt, event._tree.b_out_Lepton3_eta)
            sf_lepton3_errStat = self._worker_el_StatErr.getSF(event._tree.b_out_Lepton3_pdgId, event._tree.b_out_Lepton3_pt, event._tree.b_out_Lepton3_eta)
            sf_lepton3_errSyst = self._worker_el_SystErr.getSF(event._tree.b_out_Lepton3_pdgId, event._tree.b_out_Lepton3_pt, event._tree.b_out_Lepton3_eta)
            sf_el *= sf_lepton3
            sf_el_err2 += sf_lepton3_err**2
            sf_el_err2Stat += sf_lepton3_errStat**2
            sf_el_err2Syst += sf_lepton3_errSyst**2
        elif abs(event._tree.b_out_Lepton3_pdgId) == 13:
            sf_lepton3 = self._worker_mu.getSF(event._tree.b_out_Lepton3_pdgId, event._tree.b_out_Lepton3_pt, event._tree.b_out_Lepton3_eta)
            sf_lepton3_err = self._worker_mu.getSFErr(event._tree.b_out_Lepton3_pdgId, event._tree.b_out_Lepton3_pt, event._tree.b_out_Lepton3_eta)
            sf_lepton3_errStat = self._worker_mu_Stat.getSFErr(event._tree.b_out_Lepton3_pdgId, event._tree.b_out_Lepton3_pt, event._tree.b_out_Lepton3_eta)
            sf_lepton3_errSyst = self._worker_mu_Syst.getSFErr(event._tree.b_out_Lepton3_pdgId, event._tree.b_out_Lepton3_pt, event._tree.b_out_Lepton3_eta)
            sf_mu *= sf_lepton3
            sf_mu_err2 += sf_lepton3_err**2
            sf_mu_err2Stat += sf_lepton3_errStat**2
            sf_mu_err2Syst += sf_lepton3_errSyst**2

        sf_lepton_total = sf_lepton1 * sf_lepton2 * sf_lepton3
        sf_lepton_err_total = sqrt( sf_lepton1_err**2 + sf_lepton2_err**2 + sf_lepton3_err**2 )
        sf_lepton_errStat_total = sqrt( sf_lepton1_errStat**2 + sf_lepton2_errStat**2 + sf_lepton3_errStat**2 )
        sf_lepton_errSyst_total = sqrt( sf_lepton1_errSyst**2 + sf_lepton2_errSyst**2 + sf_lepton3_errSyst**2 )
        sf_el_err = sqrt(sf_el_err2)
        sf_el_errStat = sqrt(sf_el_err2Stat)
        sf_el_errSyst = sqrt(sf_el_err2Syst)
        sf_mu_err = sqrt(sf_mu_err2)
        sf_mu_errStat = sqrt(sf_mu_err2Stat)
        sf_mu_errSyst = sqrt(sf_mu_err2Syst)

        self.out.fillBranch("FakeLepton_SF", sf_lepton_total)
        self.out.fillBranch("FakeLepton_SFerr", sf_lepton_err_total)
        self.out.fillBranch("FakeLepton_SFerrStat", sf_lepton_errStat_total)
        self.out.fillBranch("FakeLepton_SFerrSyst", sf_lepton_errSyst_total)
        #self.out.fillBranch("FakeLepton1_SF", sf_lepton1)
        #self.out.fillBranch("FakeLepton1_SFerr", sf_lepton1_err)
        #self.out.fillBranch("FakeLepton2_SF", sf_lepton2)
        #self.out.fillBranch("FakeLepton2_SFerr", sf_lepton2_err)
        #self.out.fillBranch("FakeLepton3_SF", sf_lepton3)
        #self.out.fillBranch("FakeLepton3_SFerr", sf_lepton3_err)
        self.out.fillBranch("FakeElectron_SF", sf_el)
        self.out.fillBranch("FakeElectron_SFerr", sf_el_err)
        self.out.fillBranch("FakeElectron_SFerrStat", sf_el_errStat)
        self.out.fillBranch("FakeElectron_SFerrSyst", sf_el_errSyst)
        self.out.fillBranch("FakeMuon_SF", sf_mu)
        self.out.fillBranch("FakeMuon_SFerr", sf_mu_err)
        self.out.fillBranch("FakeMuon_SFerrStat", sf_mu_errStat)
        self.out.fillBranch("FakeMuon_SFerrSyst", sf_mu_errSyst)
        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

lepSF_forfake = lambda : lepSFProducerForFake()
