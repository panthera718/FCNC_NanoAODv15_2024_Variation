import os

import correctionlib
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class puWeightProducer(Module):
    """
    Official 2024 pileup reweighting producer using the LUM POG correctionlib
    JSON for Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15.

    The input variable is NanoAOD Pileup_nTrueInt.  The output branch names are
    kept compatible with the old workflow: puWeight, puWeightUp, puWeightDown.
    """

    def __init__(self, *args, **kwargs):
        cmssw_base = os.environ.get("CMSSW_BASE", "")
        self.era_set = kwargs.get("eraSet", "CDEFGHI")
        self.name = kwargs.get("name", "puWeight")
        self.doSysVar = bool(kwargs.get("doSysVar", True))
        self.ntrueint_branch = kwargs.get("nTrueIntBranch", "Pileup_nTrueInt")

        default_local = os.path.join(
            cmssw_base,
            "src/TZWi/TopAnalysis/data/pileup/2024/puWeights_%s.json.gz" % self.era_set,
        )
        default_cvmfs = (
            "/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/"
            "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/"
            "puWeights_%s.json.gz" % self.era_set
        )
        self.json_path = kwargs.get("jsonPath", default_local if os.path.exists(default_local) else default_cvmfs)
        self.correction_name = kwargs.get("correctionName", "Collisions24_%s_goldenJSON" % self.era_set)

        if not os.path.exists(self.json_path):
            raise RuntimeError("LUM pileup JSON file not found: %s" % self.json_path)

        print("[puWeightProducer] Loading official 2024 pileup weights: %s" % self.json_path)
        self.cset = correctionlib.CorrectionSet.from_file(self.json_path)
        self.corr = self.cset[self.correction_name]

    def beginJob(self):
        pass

    def endJob(self):
        for attr in ("corr", "cset"):
            if hasattr(self, attr):
                delattr(self, attr)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self.name, "F")
        if self.doSysVar:
            self.out.branch(self.name + "Up", "F")
            self.out.branch(self.name + "Down", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def _ntrueint(self, event):
        if hasattr(event, self.ntrueint_branch):
            return float(getattr(event, self.ntrueint_branch))
        if hasattr(event, "Pileup_nPU"):
            return float(getattr(event, "Pileup_nPU"))
        return None

    def _eval(self, ntrueint, variation):
        if ntrueint is None:
            return 1.0
        try:
            return float(self.corr.evaluate(float(ntrueint), variation))
        except Exception:
            return 1.0

    def analyze(self, event):
        ntrueint = self._ntrueint(event)
        nominal = self._eval(ntrueint, "nominal")
        self.out.fillBranch(self.name, nominal)
        if self.doSysVar:
            self.out.fillBranch(self.name + "Up", self._eval(ntrueint, "up"))
            self.out.fillBranch(self.name + "Down", self._eval(ntrueint, "down"))
        return True


puWeight_2024 = lambda: puWeightProducer(eraSet="CDEFGHI")
puWeight_NewUL2024 = lambda: puWeightProducer(eraSet="CDEFGHI")
puWeight_BCDEFGHI_2024 = lambda: puWeightProducer(eraSet="BCDEFGHI")
puWeight_CDEFGHI_2024 = lambda: puWeightProducer(eraSet="CDEFGHI")
puWeight_B_2024 = lambda: puWeightProducer(eraSet="B")
puWeight_C_2024 = lambda: puWeightProducer(eraSet="C")
puWeight_D_2024 = lambda: puWeightProducer(eraSet="D")
puWeight_E_2024 = lambda: puWeightProducer(eraSet="E")
puWeight_F_2024 = lambda: puWeightProducer(eraSet="F")
puWeight_G_2024 = lambda: puWeightProducer(eraSet="G")
puWeight_H_2024 = lambda: puWeightProducer(eraSet="H")
puWeight_I_2024 = lambda: puWeightProducer(eraSet="I")
