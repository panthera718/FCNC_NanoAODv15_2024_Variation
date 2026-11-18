import ROOT
import os
import math
import copy
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class _WeightCalculatorFromHistogramPy(object):
    """Pure-Python replacement for WeightCalculatorFromHistogram used here only
    in the simple histogram-lookup mode: WeightCalculatorFromHistogram(hist).
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
            if not tf:
                print("WARNING! File %s cannot be opened. Skipping this scale factor " % filename)
                continue

            obj = tf.Get(histname)
            if not obj:
                print("ERROR! Histogram %s not in file %s. Not considering this SF. " % (histname, filename))
                tf.Close()
                continue

            hist = obj.Clone(("eff_" + histname))
            hist.SetDirectory(0)
            print("Loading histogram %s from file %s... " % (histname, filename))
            self.effmaps_.append(hist)
            tf.Close()

    def getJV(self, eta, phi):
        out = 1.0
        x = eta
        y = phi
        for hist in self.effmaps_:
            wc = _WeightCalculatorFromHistogramPy(hist)
            out *= wc.getWeight(x, y)
        return float(out)


class JetVetoProducer(Module):

    def __init__(self):

        jveto_fN = ["hotjets-UL18.root"]
        jveto_h = ["h2hot_ul18_plus_hem1516_and_hbp2m1"]

        jveto_f = ["%s/src/TZWi/TopAnalysis/data/JetVeto/" % os.environ['CMSSW_BASE'] + f for f in jveto_fN]

        self.jveto_f = list(jveto_f)
        self.jveto_h = list(jveto_h)

    def beginJob(self):
        self._worker_jveto = _JetVetoCorrectorPy(self.jveto_f, self.jveto_h)

    def endJob(self):
        if hasattr(self, '_worker_jveto'):
            del self._worker_jveto

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
#        self.out.branch("JetVeto", "F", lenVar="nGoodJet")
        self.out.branch("JetVetoSum", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jetveto = []
        jetvetosum = 0.

        for jindex in range(event._tree.b_out_nGoodJet):
            jeta = (event._tree.b_out_GoodJet_eta)[jindex]
            jphi = (event._tree.b_out_GoodJet_phi)[jindex]
            jetvetoElement = self._worker_jveto.getJV(jeta, jphi)
            jetveto.append(jetvetoElement)
            jetvetosum += jetvetoElement

#        setattr(event._tree, "b_out_JetVeto", jetveto)
        setattr(event._tree, "b_out_JetVetoSum", jetvetosum)

#        self.out.fillBranch("JetVeto", jetveto)
        self.out.fillBranch("JetVetoSum", jetvetosum)

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

JVProducer = lambda: JetVetoProducer()
