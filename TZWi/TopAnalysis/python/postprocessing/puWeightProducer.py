from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
import os

ROOT.PyConfig.IgnoreCommandLineOptions = True


class _WeightCalculatorFromHistogramPy(object):
    """
    Pure-Python port of TZWi/TopAnalysis WeightCalculatorFromHistogram.
    Supports the same two use-cases as the C++ class:
      1) direct lookup from a histogram
      2) ratio(hist, targethist) lookup with optional normalization and
         large-weight mitigation
    """

    def __init__(self, hist=None, targethist=None, norm=True, fixLargeWeights=True, verbose=False):
        self.histogram_ = None
        self.refvals_ = []
        self.targetvals_ = []
        self.verbose_ = bool(verbose)
        self.norm_ = bool(norm)

        if hist is None:
            return

        if targethist is None:
            self.histogram_ = hist
            return

        if hist.GetNcells() != targethist.GetNcells():
            print("ERROR! Numerator and denominator histograms have different number of bins!")
            self.histogram_ = None
            return

        for i in range(int(hist.GetNcells())):
            self.refvals_.append(float(hist.GetBinContent(i)))
            self.targetvals_.append(float(targethist.GetBinContent(i)))

        self.histogram_ = self._ratio(hist, targethist, bool(fixLargeWeights))

    def _loadVals(self, hist, norm=True):
        nbins = int(hist.GetNcells())
        vals = []
        for i in range(nbins):
            bc = float(hist.GetBinContent(i))
            vals.append(max(bc, 0.0))
        if self.verbose_:
            print("Normalization of %s: %s" % (hist.GetName(), hist.Integral()))
        if norm:
            integral = float(hist.Integral())
            scale = (1.0 / integral) if integral != 0.0 else 0.0
            for i in range(nbins):
                vals[i] *= scale
        return vals

    def _checkIntegral(self, wgt1, wgt2):
        myint = 0.0
        refint = 0.0
        for i in range(len(wgt1)):
            myint += float(wgt1[i]) * float(self.refvals_[i])
            refint += float(wgt2[i]) * float(self.refvals_[i])
        if refint == 0.0:
            return 0.0
        return (myint - refint) / refint

    def _fixLargeWeights(self, weights, maxshift=0.0025, hardmax=3.0):
        if not weights:
            return
        maxw = min(max(weights), 5.0)
        cropped = []
        while maxw > hardmax:
            cropped = [min(maxw, w) for w in weights]
            shift = self._checkIntegral(cropped, weights)
            if self.verbose_:
                print("For maximum weight %s: integral relative change: %s" % (maxw, shift))
            if abs(shift) > maxshift:
                break
            maxw *= 0.95
        maxw /= 0.95
        if len(cropped) > 0:
            for i in range(len(weights)):
                cropped[i] = min(maxw, weights[i])
            normshift = self._checkIntegral(cropped, weights)
            for i in range(len(weights)):
                weights[i] = cropped[i] * (1.0 - normshift)

    def _ratio(self, hist, targethist, fixLargeWgts):
        ret = hist.Clone("hweights")
        ret.SetDirectory(0)

        vals = self._loadVals(hist, self.norm_)
        targetvals = self._loadVals(targethist, self.norm_)
        weights = []
        nbins = len(vals)
        if self.verbose_:
            print("Weights for variable %s with a number of bins equal to %d:" % (hist.GetName(), nbins))
        for i in range(nbins):
            weight = (targetvals[i] / vals[i]) if vals[i] != 0 else 1.0
            if self.verbose_:
                print("%.3g" % weight, end=" ")
            weights.append(weight)
        if self.verbose_:
            print(".")
        if fixLargeWgts:
            self._fixLargeWeights(weights)
        if self.verbose_:
            print("Final weights:")
        for i in range(len(weights)):
            ret.SetBinContent(i, float(weights[i]))
            ret.SetBinError(i, 0.0)
            if self.verbose_:
                print("%.3g" % weights[i], end=" ")
        if self.verbose_:
            print(".")
        return ret

    def getWeight(self, x, y=0):
        if self.histogram_ is None:
            print("ERROR! The weights input histogram is not loaded. Returning weight 0!")
            return 0.0
        if not self.histogram_.InheritsFrom("TH2"):
            binx = max(1, min(int(self.histogram_.GetNbinsX()), int(self.histogram_.GetXaxis().FindBin(x))))
            return float(self.histogram_.GetBinContent(binx))
        binx = max(1, min(int(self.histogram_.GetNbinsX()), int(self.histogram_.GetXaxis().FindBin(x))))
        biny = max(1, min(int(self.histogram_.GetNbinsY()), int(self.histogram_.GetYaxis().FindBin(y))))
        return float(self.histogram_.GetBinContent(binx, biny))

    def getWeightErr(self, x, y=0):
        if self.histogram_ is None:
            print("ERROR! The weights input histogram is not loaded. Returning weight error 1!")
            return 1.0
        if not self.histogram_.InheritsFrom("TH2"):
            binx = max(1, min(int(self.histogram_.GetNbinsX()), int(self.histogram_.GetXaxis().FindBin(x))))
            return float(self.histogram_.GetBinError(binx))
        binx = max(1, min(int(self.histogram_.GetNbinsX()), int(self.histogram_.GetXaxis().FindBin(x))))
        biny = max(1, min(int(self.histogram_.GetNbinsY()), int(self.histogram_.GetYaxis().FindBin(y))))
        return float(self.histogram_.GetBinError(binx, biny))


class puWeightProducer(Module):
    def __init__(self,
                 myfile,
                 targetfile_central,
                 targetfile_plus,
                 targetfile_minus,
                 myhist="pileup",
                 targethist="pileup",
                 name="puWeight",
                 norm=True,
                 verbose=False,
                 nvtx_var="Pileup_nTrueInt",
                 doSysVar=True):

        # Compatibility layer:
        # The uploaded file mixes two positional calling conventions.
        # Legacy calls pass (..., "pu_mc", "pileup", verbose=False, ...)
        # while NewUL2018 passes explicit central/up/down ROOT files. Detect the
        # legacy non-ROOT pattern and reinterpret the arguments accordingly.
        legacy_positional = (
            isinstance(targetfile_plus, str) and isinstance(targetfile_minus, str)
            and not targetfile_plus.lower().endswith(".root")
            and not targetfile_minus.lower().endswith(".root")
        )
        if legacy_positional:
            myhist = targetfile_plus
            targethist = targetfile_minus
            targetfile_plus = None
            targetfile_minus = None
            if bool(doSysVar):
                if verbose:
                    print("[puWeightProducer] Legacy constructor pattern detected; using central histogram for Up/Down as well.")
            # Preserve branch presence when caller requested doSysVar=True.

        self.targeth = self.loadHisto(targetfile_central, targethist)
        self.doSysVar = bool(doSysVar)
        if self.doSysVar:
            if targetfile_plus is not None and targetfile_minus is not None:
                self.targeth_plus = self.loadHisto(targetfile_plus, targethist)
                self.targeth_minus = self.loadHisto(targetfile_minus, targethist)
            else:
                self.targeth_plus = self.targeth.Clone(self.targeth.GetName() + "__plus")
                self.targeth_plus.SetDirectory(0)
                self.targeth_minus = self.targeth.Clone(self.targeth.GetName() + "__minus")
                self.targeth_minus.SetDirectory(0)

        self.fixLargeWeights = True  # temporary fix
        if myfile != "auto":
            self.autoPU = False
            self.myh = self.loadHisto(myfile, myhist)
        else:
            self.fixLargeWeights = False  # AR: it seems to crash with it, to be debugged
            self.autoPU = True
            ROOT.gROOT.cd()
            self.myh = self.targeth.Clone("autoPU")
            self.myh.Reset()

        self.name = name
        self.norm = bool(norm)
        self.verbose = bool(verbose)
        self.nvtxVar = nvtx_var

    def loadHisto(self, filename, hname):
        tf = ROOT.TFile.Open(filename)
        if not tf or tf.IsZombie():
            raise RuntimeError("Cannot open ROOT file: %s" % filename)

        hist = tf.Get(hname)
        if not hist:
            tf.Close()
            raise RuntimeError("Histogram '%s' not found in file: %s" % (hname, filename))

        hist = hist.Clone("%s__%s" % (hname, os.path.basename(filename).replace(".", "_")))
        hist.SetDirectory(0)
        tf.Close()
        return hist

    def beginJob(self):
        pass

    def endJob(self):
        if hasattr(self, 'targeth'):
            del self.targeth
        if hasattr(self, 'targeth_plus'):
            del self.targeth_plus
        if hasattr(self, 'targeth_minus'):
            del self.targeth_minus
        if hasattr(self, 'myh'):
            del self.myh

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        if self.autoPU:
            self.myh.Reset()
            print("Computing PU profile for this file")
            ROOT.gROOT.cd()
            inputFile.Get("Events").Project("autoPU", self.nvtxVar)
            if outputFile:
                outputFile.cd()
                self.myh.Write()

        self._worker = _WeightCalculatorFromHistogramPy(
            self.myh, self.targeth, self.norm, self.fixLargeWeights, self.verbose
        )
        self.out = wrappedOutputTree
        self.out.branch(self.name, "F")
        if self.doSysVar:
            self._worker_plus = _WeightCalculatorFromHistogramPy(
                self.myh, self.targeth_plus, self.norm, self.fixLargeWeights, self.verbose
            )
            self._worker_minus = _WeightCalculatorFromHistogramPy(
                self.myh, self.targeth_minus, self.norm, self.fixLargeWeights, self.verbose
            )
            self.out.branch(self.name + "Up", "F")
            self.out.branch(self.name + "Down", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        if hasattr(self, '_worker'):
            del self._worker
        if hasattr(self, '_worker_plus'):
            del self._worker_plus
        if hasattr(self, '_worker_minus'):
            del self._worker_minus

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        weight = 1.0
        weight_plus = 1.0
        weight_minus = 1.0
        if hasattr(event, self.nvtxVar):
            nvtx = int(getattr(event, self.nvtxVar))
            weight = self._worker.getWeight(nvtx) if nvtx < self.myh.GetNbinsX() else 1.0
            if self.doSysVar:
                weight_plus = self._worker_plus.getWeight(nvtx) if nvtx < self.myh.GetNbinsX() else 1.0
                weight_minus = self._worker_minus.getWeight(nvtx) if nvtx < self.myh.GetNbinsX() else 1.0
        self.out.fillBranch(self.name, weight)
        if self.doSysVar:
            self.out.fillBranch(self.name + "Up", weight_plus)
            self.out.fillBranch(self.name + "Down", weight_minus)
        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

# NewUL2018
pufile_dataNewUL2018_central = "%s/src/TZWi/TopAnalysis/data/pileup/PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins_New.root" % os.environ['CMSSW_BASE']
pufile_dataNewUL2018_plus = "%s/src/TZWi/TopAnalysis/data/pileup/PileupHistogram-goldenJSON-13tev-2018-72400ub-99bins_New.root" % os.environ['CMSSW_BASE']
pufile_dataNewUL2018_minus = "%s/src/TZWi/TopAnalysis/data/pileup/PileupHistogram-goldenJSON-13tev-2018-66000ub-99bins_New.root" % os.environ['CMSSW_BASE']
pufile_mcNewUL2018 = "%s/src/TZWi/TopAnalysis/data/pileup/mcPileupUL2018.root" % os.environ['CMSSW_BASE']
puWeight_NewUL2018 = lambda: puWeightProducer(pufile_mcNewUL2018,
                                           pufile_dataNewUL2018_central,
                                           pufile_dataNewUL2018_plus,
                                           pufile_dataNewUL2018_minus,
                                           "pu_mc",
                                           "pileup_new",
                                           verbose=False,
                                           doSysVar=True)
