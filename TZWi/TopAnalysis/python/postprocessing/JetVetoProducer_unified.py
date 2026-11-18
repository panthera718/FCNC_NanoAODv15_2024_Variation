import ROOT
import os
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection  # kept for compatibility
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class _WeightCalculatorFromHistogramPy(object):
    """Pure-Python replacement for WeightCalculatorFromHistogram in simple
    histogram-lookup mode.
    """

    def __init__(self, histogram, verbose=False):
        self.histogram_ = histogram
        self.verbose_ = verbose

    def getWeight(self, x, y=0):
        if self.histogram_ is None:
            print("ERROR! The weights input histogram is not loaded. Returning weight 0!")
            return 0.0

        if not self.histogram_.InheritsFrom("TH2"):
            binx = max(1, min(self.histogram_.GetNbinsX(), self.histogram_.GetXaxis().FindBin(x)))
            return float(self.histogram_.GetBinContent(binx))

        binx = max(1, min(self.histogram_.GetNbinsX(), self.histogram_.GetXaxis().FindBin(x)))
        biny = max(1, min(self.histogram_.GetNbinsY(), self.histogram_.GetYaxis().FindBin(y)))
        return float(self.histogram_.GetBinContent(binx, biny))


class _JetVetoCorrectorPy(object):
    def __init__(self, files, histos):
        self.effmaps_ = []
        if len(files) != len(histos):
            print("ERROR! There should be one histogram per input file! Returning 0 as SF.")
            return

        for filename, histname in zip(files, histos):
            tf = ROOT.TFile.Open(filename, "read")
            if not tf or tf.IsZombie():
                print("WARNING! File %s cannot be opened. Skipping this scale factor " % filename)
                continue

            obj = tf.Get(histname)
            if not obj:
                print("ERROR! Histogram %s not in file %s. Not considering this SF. " % (histname, filename))
                tf.Close()
                continue

            hist = obj.Clone(("eff_" + histname).replace("/", "_"))
            hist.SetDirectory(0)
            print("Loading histogram %s from file %s... " % (histname, filename))
            self.effmaps_.append(hist)
            tf.Close()

    def getJV(self, eta, phi):
        out = 1.0
        for hist in self.effmaps_:
            wc = _WeightCalculatorFromHistogramPy(hist)
            out *= wc.getWeight(eta, phi)
        return float(out)


class JetVetoProducer(Module, object):
    """
    Unified pure-Python JetVeto producer.

    Modes:
      - central: nominal branches
      - up:      jetVar Up branches
      - down:    jetVar Down branches

    All three modes use the same 2018 hot-jet map:
      hotjets-UL18.root / h2hot_ul18_plus_hem1516_and_hbp2m1
    """

    VALID_DIRECTIONS = {"central", "up", "down"}

    def __init__(self, *args, **kwargs):
        self.jetVar = kwargs.get("jetVar")
        self.direction = kwargs.get("direction", "central")
        if self.direction not in self.VALID_DIRECTIONS:
            raise ValueError("direction must be one of %s" % sorted(self.VALID_DIRECTIONS))

        base_dir = "%s/src/TZWi/TopAnalysis/data/JetVeto/" % os.environ["CMSSW_BASE"]

        if self.direction != "central" and not self.jetVar:
            raise ValueError("jetVar is required for direction='%s'" % self.direction)

        jveto_files = [base_dir + "hotjets-UL18.root"]
        self.jveto_f = list(jveto_files)
        self.jveto_h = ["h2hot_ul18_plus_hem1516_and_hbp2m1"]

    def beginJob(self):
        self._worker_jveto = _JetVetoCorrectorPy(self.jveto_f, self.jveto_h)

    def endJob(self):
        if hasattr(self, "_worker_jveto"):
            del self._worker_jveto

    def _branch_name(self):
        if self.direction == "central":
            return "JetVetoSum"
        suffix = "Up" if self.direction == "up" else "Down"
        return "JetVetoSum_%s%s" % (self.jetVar, suffix)

    def _ngoodjet_attr(self):
        if self.direction == "central":
            return "b_out_nGoodJet"
        suffix = "Up" if self.direction == "up" else "Down"
        return "b_out_nGoodJet_%s%s" % (self.jetVar, suffix)

    def _eta_attr(self):
        if self.direction == "central":
            return "b_out_GoodJet_eta"
        suffix = "Up" if self.direction == "up" else "Down"
        return "b_out_GoodJet_eta_%s%s" % (self.jetVar, suffix)

    def _phi_attr(self):
        if self.direction == "central":
            return "b_out_GoodJet_phi"
        suffix = "Up" if self.direction == "up" else "Down"
        return "b_out_GoodJet_phi_%s%s" % (self.jetVar, suffix)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self._branch_name(), "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        tree = event._tree
        jetveto = []
        jetvetosum = 0.0

        n_attr = self._ngoodjet_attr()
        eta_attr = self._eta_attr()
        phi_attr = self._phi_attr()
        out_attr = self._branch_name()

        n_goodjet = int(getattr(tree, n_attr))
        jet_etas = getattr(tree, eta_attr)
        jet_phis = getattr(tree, phi_attr)

        for jindex in range(n_goodjet):
            jeta = jet_etas[jindex]
            jphi = jet_phis[jindex]
            jetveto_element = self._worker_jveto.getJV(jeta, jphi)
            jetveto.append(jetveto_element)
            jetvetosum += jetveto_element

        setattr(tree, "b_out_%s" % out_attr, jetvetosum)
        self.out.fillBranch(out_attr, jetvetosum)
        return True


# Central lambda
JVProducer = lambda: JetVetoProducer(direction="central")


# Unified-file lambdas for Up/Down usage with the same imported module.
# Example shell usage after unification:
#   -I TZWi.TopAnalysis.postprocessing.JetVetoProducer JVProducer
#   -I TZWi.TopAnalysis.postprocessing.JetVetoProducer JVProducerUp_${jetVar}
#   -I TZWi.TopAnalysis.postprocessing.JetVetoProducer JVProducerDown_${jetVar}

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
    return lambda direction=direction, jetVar=jetVar: JetVetoProducer(direction=direction, jetVar=jetVar)


for _jetVar in _JETVARS:
    globals()["JVProducerUp_%s" % _jetVar] = _make_lambda("up", _jetVar)
    globals()["JVProducerDown_%s" % _jetVar] = _make_lambda("down", _jetVar)
