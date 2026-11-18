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
    It reproduces the behavior relevant for lepRecoSFProducer:
      - load TH2 histograms from ROOT files
      - clip pt to the same boundaries as the C++ code
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

    def getSFErr(self, pdgid, pt, eta):
        out = 0.0
        pt = self._clip_pt(pdgid, pt)
        x = abs(float(eta)) if abs(pdgid) == 13 else float(eta)
        y = pt
        for hist in self.effmaps_:
            wc = _HistogramWeightCalculator(hist)
            out += wc.getWeightErr(x, y) ** 2.0
        return float(sqrt(out))

    def run(self):
        self.ret_ = []
        for iL in range(self.nLep_):
            self.ret_.append(self.getSF(self.Lep_pdgId_[iL], self.Lep_pt_[iL], self.Lep_eta_[iL]))
        return self.ret_


class lepRecoSFProducer(Module):
    def __init__(self):
        el_below_f = ["egammaEffi_ptBelow20.txt_EGM2D_UL2018.root"]
        el_above_f = ["egammaEffi_ptAbove20.txt_EGM2D_UL2018.root"]
        el_below_h = ["EGamma_SF2D"]
        el_above_h = ["EGamma_SF2D"]

        base = "%s/src/TZWi/TopAnalysis/data/leptonRecoSF/" % os.environ['CMSSW_BASE']
        self.el_below_f = [base + f for f in el_below_f]
        self.el_above_f = [base + f for f in el_above_f]
        self.el_below_h = list(el_below_h)
        self.el_above_h = list(el_above_h)

    def beginJob(self):
        self._worker_el_below = _PyLeptonEfficiencyCorrector(self.el_below_f, self.el_below_h)
        self._worker_el_above = _PyLeptonEfficiencyCorrector(self.el_above_f, self.el_above_h)

    def endJob(self):
        if hasattr(self, '_worker_el_below'):
            del self._worker_el_below
        if hasattr(self, '_worker_el_above'):
            del self._worker_el_above

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Lepton_RecoSF", "F")
        self.out.branch("Lepton_RecoSFerr", "F")
        self.out.branch("Electron_RecoSF", "F")
        self.out.branch("Electron_RecoSFerr", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def _one_electron_sf(self, pdgid, pt, eta):
        if abs(int(pdgid)) != 11:
            return 1.0, 0.0
        if pt < 20.0:
            return (
                self._worker_el_below.getSF(pdgid, pt, eta),
                self._worker_el_below.getSFErr(pdgid, pt, eta),
            )
        return (
            self._worker_el_above.getSF(pdgid, pt, eta),
            self._worker_el_above.getSFErr(pdgid, pt, eta),
        )

    def analyze(self, event):
        sf_lepton_total = 1.0
        sf_lepton_err2 = 0.0
        sf_el = 1.0
        sf_el_err2 = 0.0

        for idx in (1, 2, 3):
            pdgid = getattr(event._tree, "b_out_Lepton%d_pdgId" % idx)
            pt = getattr(event._tree, "b_out_Lepton%d_pt" % idx)
            eta = getattr(event._tree, "b_out_Lepton%d_eta" % idx)

            sf, sf_err = self._one_electron_sf(pdgid, pt, eta)
            sf_lepton_total *= sf
            sf_lepton_err2 += sf_err ** 2

            if abs(int(pdgid)) == 11:
                sf_el *= sf
                sf_el_err2 += sf_err ** 2

        self.out.fillBranch("Lepton_RecoSF", sf_lepton_total)
        self.out.fillBranch("Lepton_RecoSFerr", sqrt(sf_lepton_err2))
        self.out.fillBranch("Electron_RecoSF", sf_el)
        self.out.fillBranch("Electron_RecoSFerr", sqrt(sf_el_err2))
        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having
# them loaded when not needed
lepRecoSF = lambda: lepRecoSFProducer()
