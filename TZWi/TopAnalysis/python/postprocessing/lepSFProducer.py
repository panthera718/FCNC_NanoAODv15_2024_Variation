import os
from math import sqrt

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class _HistogramWeightCalculator(object):
    """
    Pure-Python equivalent of WeightCalculatorFromHistogram for the only
    functionality used here: reading a value or its bin error from a TH1/TH2
    with under/overflow clamped to the edge bins.
    """

    def __init__(self, histogram):
        self.histogram = histogram

    def _clamp_xbin(self, x):
        axis = self.histogram.GetXaxis()
        return max(1, min(self.histogram.GetNbinsX(), axis.FindBin(x)))

    def _clamp_ybin(self, y):
        axis = self.histogram.GetYaxis()
        return max(1, min(self.histogram.GetNbinsY(), axis.FindBin(y)))

    def getWeight(self, x, y=0.0):
        if self.histogram is None:
            print("ERROR! The weights input histogram is not loaded. Returning weight 0!")
            return 0.0
        if not self.histogram.InheritsFrom("TH2"):
            return float(self.histogram.GetBinContent(self._clamp_xbin(x)))
        return float(self.histogram.GetBinContent(self._clamp_xbin(x), self._clamp_ybin(y)))

    def getWeightErr(self, x, y=0.0):
        if self.histogram is None:
            print("ERROR! The weights input histogram is not loaded. Returning weight error 1!")
            return 1.0
        if not self.histogram.InheritsFrom("TH2"):
            return float(self.histogram.GetBinError(self._clamp_xbin(x)))
        return float(self.histogram.GetBinError(self._clamp_xbin(x), self._clamp_ybin(y)))


class _PyLeptonEfficiencyCorrector(object):
    """
    Pure-Python replacement for LeptonEfficiencyCorrector.
    It reproduces the behavior of the C++ implementation used by lepSFProducer:
      - load TH2 histograms from ROOT files
      - clip pt to the same boundaries
      - use |eta| for muons and signed eta for electrons
      - multiply SFs across all loaded maps
      - combine uncertainties in quadrature using histogram bin errors
    """

    def __init__(self, files, histos):
        self.effmaps_ = []
        self.ret_ = []
        self.nLep_ = 0
        self.Lep_pdgId_ = None
        self.Lep_pt_ = None
        self.Lep_eta_ = None

        if len(files) != len(histos):
            print("ERROR! There should be one histogram per input file! Returning 0 as SF.")
            return

        for i in range(len(files)):
            path = files[i]
            hist_name = histos[i]
            root_file = ROOT.TFile.Open(path, "READ")
            if not root_file or root_file.IsZombie():
                print("WARNING! File {0} cannot be opened. Skipping this scale factor".format(path))
                continue

            obj = root_file.Get(hist_name)
            if not obj:
                print("ERROR! Histogram {0} not in file {1}. Not considering this SF.".format(hist_name, path))
                root_file.Close()
                continue

            hist = obj.Clone("eff_{0}_{1}".format(i, hist_name.replace('/', '_')))
            hist.SetDirectory(0)
            print("Loading histogram {0} from file {1}...".format(hist_name, path))
            self.effmaps_.append(hist)
            root_file.Close()

    def setLeptons(self, nLep, lepPdgId, lepPt, lepEta):
        self.nLep_ = nLep
        self.Lep_pdgId_ = lepPdgId
        self.Lep_pt_ = lepPt
        self.Lep_eta_ = lepEta

    def _clip_pt(self, pdgid, pt):
        if abs(pdgid) == 13:
            if pt >= 120.0:
                return 119.0
            if pt <= 10.0:
                return 11.0
        else:
            if pt >= 200.0:
                return 199.0
            if pt <= 10.0:
                return 11.0
        return float(pt)

    def getSF(self, pdgid, pt, eta):
        out = 1.0
        pt = self._clip_pt(pdgid, pt)
        x = abs(float(eta)) if abs(pdgid) == 13 else float(eta)
        y = pt
        for hist in self.effmaps_:
            wc = _HistogramWeightCalculator(hist)
            out *= wc.getWeight(x, y)
        return float(out)

    def getTrigSF(self, lep1_pt, lep2_pt):
        out = 1.0
        if lep1_pt >= 500.0:
            lep1_pt = 499.0
        elif lep1_pt <= 20.0:
            lep1_pt = 19.0
        if lep2_pt >= 500.0:
            lep2_pt = 499.0
        elif lep2_pt <= 20.0:
            lep2_pt = 19.0
        x = float(lep1_pt)
        y = float(lep2_pt)
        for hist in self.effmaps_:
            wc = _HistogramWeightCalculator(hist)
            out *= wc.getWeight(x, y)
        return float(out)

    def getSFErr(self, pdgid, pt, eta):
        out = 0.0
        pt = self._clip_pt(pdgid, pt)
        x = abs(float(eta)) if abs(pdgid) == 13 else float(eta)
        y = pt
        for hist in self.effmaps_:
            wc = _HistogramWeightCalculator(hist)
            out += wc.getWeightErr(x, y) ** 2.0
        return float(sqrt(out))

    def getTrigSFErr(self, lep1_pt, lep2_pt):
        out = 0.0
        if lep1_pt >= 500.0:
            lep1_pt = 499.0
        elif lep1_pt <= 20.0:
            lep1_pt = 19.0
        if lep2_pt >= 500.0:
            lep2_pt = 499.0
        elif lep2_pt <= 20.0:
            lep2_pt = 19.0
        x = float(lep1_pt)
        y = float(lep2_pt)
        for hist in self.effmaps_:
            wc = _HistogramWeightCalculator(hist)
            out += wc.getWeightErr(x, y) ** 2.0
        return float(sqrt(out))

    def run(self):
        self.ret_ = []
        for iL in range(self.nLep_):
            self.ret_.append(self.getSF(self.Lep_pdgId_[iL], self.Lep_pt_[iL], self.Lep_eta_[iL]))
        return self.ret_


class lepSFProducer(Module):
    def __init__(self):
        mu_f = ["NUM_LeptonMvaMedium_DEN_TrackerMuons_abseta_pt.root"]
        mu_h = ["NUM_LeptonMvaMedium_DEN_TrackerMuons_abseta_pt"]
        mu_h_stat = ["NUM_LeptonMvaMedium_DEN_TrackerMuons_abseta_pt_stat"]
        mu_h_syst = ["NUM_LeptonMvaMedium_DEN_TrackerMuons_abseta_pt_syst"]
        el_f = ["egammaEffi.txt_EGM2D.root"]
        el_h = ["EGamma_SF2D"]
        el_h_statErr = ["stat"]
        el_h_systErr = ["sys"]

        base = "%s/src/TZWi/TopAnalysis/data/leptonSF/" % os.environ['CMSSW_BASE']
        self.mu_f = [base + f for f in mu_f]
        self.mu_h = list(mu_h)
        self.mu_h_stat = list(mu_h_stat)
        self.mu_h_syst = list(mu_h_syst)
        self.el_f = [base + f for f in el_f]
        self.el_h = list(el_h)
        self.el_h_statErr = list(el_h_statErr)
        self.el_h_systErr = list(el_h_systErr)

    def beginJob(self):
        self._worker_el = _PyLeptonEfficiencyCorrector(self.el_f, self.el_h)
        self._worker_el_StatErr = _PyLeptonEfficiencyCorrector(self.el_f, self.el_h_statErr)
        self._worker_el_SystErr = _PyLeptonEfficiencyCorrector(self.el_f, self.el_h_systErr)
        self._worker_mu = _PyLeptonEfficiencyCorrector(self.mu_f, self.mu_h)
        self._worker_mu_Stat = _PyLeptonEfficiencyCorrector(self.mu_f, self.mu_h_stat)
        self._worker_mu_Syst = _PyLeptonEfficiencyCorrector(self.mu_f, self.mu_h_syst)

    def endJob(self):
        # Explicit cleanup keeps ownership simple for PyROOT-held histograms.
        if hasattr(self, '_worker_el'):
            del self._worker_el
        if hasattr(self, '_worker_el_StatErr'):
            del self._worker_el_StatErr
        if hasattr(self, '_worker_el_SystErr'):
            del self._worker_el_SystErr
        if hasattr(self, '_worker_mu'):
            del self._worker_mu
        if hasattr(self, '_worker_mu_Stat'):
            del self._worker_mu_Stat
        if hasattr(self, '_worker_mu_Syst'):
            del self._worker_mu_Syst

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Lepton_SF", "F")
        self.out.branch("Lepton_SFerr", "F")
        self.out.branch("Lepton_SFerrStat", "F")
        self.out.branch("Lepton_SFerrSyst", "F")
        self.out.branch("Electron_SF", "F")
        self.out.branch("Electron_SFerr", "F")
        self.out.branch("Electron_SFerrStat", "F")
        self.out.branch("Electron_SFerrSyst", "F")
        self.out.branch("Muon_SF", "F")
        self.out.branch("Muon_SFerr", "F")
        self.out.branch("Muon_SFerrStat", "F")
        self.out.branch("Muon_SFerrSyst", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def _process_one_lepton(self, pdgid, pt, eta):
        sf = 1.0
        sf_err = 0.0
        sf_err_stat = 0.0
        sf_err_syst = 0.0
        flavor = None

        apdgid = abs(int(pdgid))
        if apdgid == 11:
            sf = self._worker_el.getSF(pdgid, pt, eta)
            sf_err = self._worker_el.getSFErr(pdgid, pt, eta)
            sf_err_stat = self._worker_el_StatErr.getSF(pdgid, pt, eta)
            sf_err_syst = self._worker_el_SystErr.getSF(pdgid, pt, eta)
            flavor = "electron"
        elif apdgid == 13:
            sf = self._worker_mu.getSF(pdgid, pt, eta)
            sf_err = self._worker_mu.getSFErr(pdgid, pt, eta)
            sf_err_stat = self._worker_mu_Stat.getSFErr(pdgid, pt, eta)
            sf_err_syst = self._worker_mu_Syst.getSFErr(pdgid, pt, eta)
            flavor = "muon"

        return flavor, sf, sf_err, sf_err_stat, sf_err_syst

    def analyze(self, event):
        sf_lepton_total = 1.0
        sf_lepton_err2 = 0.0
        sf_lepton_err2_stat = 0.0
        sf_lepton_err2_syst = 0.0

        sf_el = 1.0
        sf_el_err2 = 0.0
        sf_el_err2_stat = 0.0
        sf_el_err2_syst = 0.0

        sf_mu = 1.0
        sf_mu_err2 = 0.0
        sf_mu_err2_stat = 0.0
        sf_mu_err2_syst = 0.0

        for idx in (1, 2, 3):
            pdgid = getattr(event._tree, "b_out_Lepton%d_pdgId" % idx)
            pt = getattr(event._tree, "b_out_Lepton%d_pt" % idx)
            eta = getattr(event._tree, "b_out_Lepton%d_eta" % idx)

            flavor, sf, sf_err, sf_err_stat, sf_err_syst = self._process_one_lepton(pdgid, pt, eta)

            sf_lepton_total *= sf
            sf_lepton_err2 += sf_err ** 2
            sf_lepton_err2_stat += sf_err_stat ** 2
            sf_lepton_err2_syst += sf_err_syst ** 2

            if flavor == "electron":
                sf_el *= sf
                sf_el_err2 += sf_err ** 2
                sf_el_err2_stat += sf_err_stat ** 2
                sf_el_err2_syst += sf_err_syst ** 2
            elif flavor == "muon":
                sf_mu *= sf
                sf_mu_err2 += sf_err ** 2
                sf_mu_err2_stat += sf_err_stat ** 2
                sf_mu_err2_syst += sf_err_syst ** 2

        self.out.fillBranch("Lepton_SF", sf_lepton_total)
        self.out.fillBranch("Lepton_SFerr", sqrt(sf_lepton_err2))
        self.out.fillBranch("Lepton_SFerrStat", sqrt(sf_lepton_err2_stat))
        self.out.fillBranch("Lepton_SFerrSyst", sqrt(sf_lepton_err2_syst))

        self.out.fillBranch("Electron_SF", sf_el)
        self.out.fillBranch("Electron_SFerr", sqrt(sf_el_err2))
        self.out.fillBranch("Electron_SFerrStat", sqrt(sf_el_err2_stat))
        self.out.fillBranch("Electron_SFerrSyst", sqrt(sf_el_err2_syst))

        self.out.fillBranch("Muon_SF", sf_mu)
        self.out.fillBranch("Muon_SFerr", sqrt(sf_mu_err2))
        self.out.fillBranch("Muon_SFerrStat", sqrt(sf_mu_err2_stat))
        self.out.fillBranch("Muon_SFerrSyst", sqrt(sf_mu_err2_syst))
        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having
# them loaded when not needed
lepSF = lambda: lepSFProducer()
