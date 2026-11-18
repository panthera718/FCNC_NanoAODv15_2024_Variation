import math
import os

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection  # kept for compatibility
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class _HistogramWeightCalculator(object):
    """
    Pure-Python equivalent of the direct-lookup mode of
    WeightCalculatorFromHistogram used by JetPUIDSFCorrector:
      - read TH1/TH2 bin content
      - clamp under/overflow to edge bins
    """

    def __init__(self, histogram):
        self.histogram = histogram

    def getWeight(self, x, y=0.0):
        if self.histogram is None:
            print("ERROR! The weights input histogram is not loaded. Returning weight 0!")
            return 0.0

        if not self.histogram.InheritsFrom("TH2"):
            binx = max(1, min(int(self.histogram.GetNbinsX()), int(self.histogram.GetXaxis().FindBin(x))))
            return float(self.histogram.GetBinContent(binx))

        binx = max(1, min(int(self.histogram.GetNbinsX()), int(self.histogram.GetXaxis().FindBin(x))))
        biny = max(1, min(int(self.histogram.GetNbinsY()), int(self.histogram.GetYaxis().FindBin(y))))
        return float(self.histogram.GetBinContent(binx, biny))


class _PyJetPUIDSFCorrector(object):
    """
    Pure-Python replacement for JetPUIDSFCorrector.
    It reproduces the C++ behavior:
      - load TH2 histograms from ROOT files
      - for each map, evaluate getWeight(pt, eta)
      - multiply all loaded-map contributions together
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
        # Kept only for interface compatibility with the original C++ class.
        self.nLep_ = nLep
        self.Lep_pdgId_ = lepPdgId
        self.Lep_pt_ = lepPt
        self.Lep_eta_ = lepEta

    def getJPUIDSF(self, pt, eta):
        out = 1.0
        x = float(pt)
        y = float(eta)
        for hist in self.effmaps_:
            wc = _HistogramWeightCalculator(hist)
            out *= wc.getWeight(x, y)
        return float(out)


class JetPUIDSFProducer(Module, object):
    """
    Unified pure-Python version of:
      - JetPUIDSFProducer.py
      - JetPUIDSFProducer_jetVarUp.py
      - JetPUIDSFProducer_jetVarDown.py

    Important compatibility behavior preserved exactly:
      - central branches are written without suffix
      - Up/Down branches are written with _<jetVar>Up / _<jetVar>Down suffix
      - Up/Down still read the *central* GoodJet branches, exactly like the
        uploaded original jetVarUp/jetVarDown implementations
    """

    VALID_DIRECTIONS = {"central", "up", "down"}

    def __init__(self, *args, **kwargs):
        self.jetVar = kwargs.get("jetVar")
        self.direction = kwargs.get("direction", "central")
        if self.direction not in self.VALID_DIRECTIONS:
            raise ValueError("direction must be one of %s" % sorted(self.VALID_DIRECTIONS))
        if self.direction in ("up", "down") and not self.jetVar:
            raise ValueError("jetVar is required for direction='%s'" % self.direction)

        jpuidsf_fN = ["PUID_106XTraining_ULRun2_EffSFandUncties_v1.root"]
        jpuidsf_h = ["h2_eff_sfUL2018_L"]
        jpuidsf_herr = ["h2_eff_sfUL2018_L_Systuncty"]

        base = "%s/src/TZWi/TopAnalysis/data/JetPUIDSF/" % os.environ['CMSSW_BASE']
        self.jpuidsf_f = [base + f for f in jpuidsf_fN]
        self.jpuidsf_h = list(jpuidsf_h)
        self.jpuidsf_herr = list(jpuidsf_herr)

    def beginJob(self):
        self._worker_jpuidsf = _PyJetPUIDSFCorrector(self.jpuidsf_f, self.jpuidsf_h)
        self._worker_jpuidsferr = _PyJetPUIDSFCorrector(self.jpuidsf_f, self.jpuidsf_herr)

    def endJob(self):
        if hasattr(self, '_worker_jpuidsf'):
            del self._worker_jpuidsf
        if hasattr(self, '_worker_jpuidsferr'):
            del self._worker_jpuidsferr

    def _branch_name(self, base):
        if self.direction == "central":
            return base
        suffix = "Up" if self.direction == "up" else "Down"
        return "%s_%s%s" % (base, self.jetVar, suffix)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self._branch_name("JetPUID_SF"), "F")
        self.out.branch(self._branch_name("JetPUID_SFerr"), "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        tree = event._tree

        jetpuid_eachsf = []
        jetpuid_eachsferr = []
        jetpuid_sf = 1.0
        jetpuid_sferr_sqsum = 0.0

        # Intentionally read the central GoodJet branches even for Up/Down,
        # matching the uploaded original JetPUIDSFProducer_jetVarUp/Down.py.
        for jetindex in range(tree.b_out_nGoodJet):
            jpt = tree.b_out_GoodJet_pt[jetindex]
            jeta = tree.b_out_GoodJet_eta[jetindex]
            if jpt < 50.0:
                jetpuidsfElement = self._worker_jpuidsf.getJPUIDSF(jpt, jeta)
                jetpuidsferrElement = self._worker_jpuidsferr.getJPUIDSF(jpt, jeta)
                jetpuid_eachsf.append(jetpuidsfElement)
                jetpuid_eachsferr.append(jetpuidsferrElement)
                jetpuid_sf *= jetpuidsfElement
                jetpuid_sferr_sqsum += jetpuidsferrElement ** 2

        jetpuid_sferr = math.sqrt(jetpuid_sferr_sqsum)

        setattr(tree, "b_out_%s" % self._branch_name("JetPUID_SF"), jetpuid_sf)
        setattr(tree, "b_out_%s" % self._branch_name("JetPUID_SFerr"), jetpuid_sferr)

        self.out.fillBranch(self._branch_name("JetPUID_SF"), jetpuid_sf)
        self.out.fillBranch(self._branch_name("JetPUID_SFerr"), jetpuid_sferr)
        return True


JPUIDSFProducer = lambda: JetPUIDSFProducer(direction="central")


_JETVARS = [
    "AbsoluteStat",
    "AbsoluteScale",
    "AbsoluteMPFBias",
    "Fragmentation",
    "SinglePionECAL",
    "SinglePionHCAL",
    "FlavorQCD",
    "TimePtEta",
    "RelativeJEREC1",
    "RelativeJEREC2",
    "RelativeJERHF",
    "RelativePtBB",
    "RelativePtEC1",
    "RelativePtEC2",
    "RelativePtHF",
    "RelativeBal",
    "RelativeSample",
    "RelativeFSR",
    "RelativeStatFSR",
    "RelativeStatEC",
    "RelativeStatHF",
    "PileUpDataMC",
    "PileUpPtRef",
    "PileUpPtBB",
    "PileUpPtEC1",
    "PileUpPtEC2",
    "PileUpPtHF",
]


def _make_lambda(direction, jetVar):
    return lambda direction=direction, jetVar=jetVar: JetPUIDSFProducer(direction=direction, jetVar=jetVar)


for _jetVar in _JETVARS:
    globals()["JPUIDSFProducerUp_%s" % _jetVar] = _make_lambda("up", _jetVar)
    globals()["JPUIDSFProducerDown_%s" % _jetVar] = _make_lambda("down", _jetVar)
