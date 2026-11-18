import os
import correctionlib

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class JetIDProducer(Module, object):
    """
    Official Run-3 AK4PUPPI jet-ID producer for 2024 NanoAODv15.

    Historical note: the old file named JetIDSFProducer_unified.py contained a
    Run-2 JetPUID scale-factor producer.  That correction should not be applied
    to the 2024 AK4PUPPI object definition used here.  This replacement only
    evaluates the official JME jet-ID decision from jetid.json.gz and writes
    validation branches.  The FCNCTriLepton producer also evaluates the same
    AK4PUPPI_Tight criterion internally for the GoodJet selection.
    """

    def __init__(self, *args, **kwargs):
        cmssw_base = os.environ.get("CMSSW_BASE", "")
        self.json_path = kwargs.get(
            "jsonPath",
            os.path.join(cmssw_base, "src/TZWi/TopAnalysis/data/JetID/2024/jetid.json.gz"),
        )
        self.correction_name = kwargs.get("correctionName", "AK4PUPPI_Tight")
        if not os.path.exists(self.json_path):
            raise RuntimeError("JME jet-ID JSON file not found: %s" % self.json_path)

        print("[JetIDProducer] Loading official 2024 jet ID: %s" % self.json_path)
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
        self.out.branch("JetID_Tight", "O", lenVar="nJet")
        self.out.branch("nJetIDTight", "i")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    @staticmethod
    def _to_int(value):
        try:
            return int(value)
        except Exception:
            try:
                return int(ord(value))
            except Exception:
                return 0

    def _passes_jet_id(self, jet):
        try:
            ch_mult = self._to_int(getattr(jet, "chMultiplicity", 0))
            ne_mult = self._to_int(getattr(jet, "neMultiplicity", 0))
            multiplicity = ch_mult + ne_mult
            return bool(
                float(
                    self.corr.evaluate(
                        float(jet.eta),
                        float(getattr(jet, "chHEF", 0.0)),
                        float(getattr(jet, "neHEF", 0.0)),
                        float(getattr(jet, "chEmEF", 0.0)),
                        float(getattr(jet, "neEmEF", 0.0)),
                        float(getattr(jet, "muEF", 0.0)),
                        int(ch_mult),
                        int(ne_mult),
                        int(multiplicity),
                    )
                )
                > 0.5
            )
        except Exception:
            return False

    def analyze(self, event):
        jetid = []
        try:
            jets = Collection(event, "Jet")
            for jet in jets:
                jetid.append(self._passes_jet_id(jet))
        except Exception:
            jetid = []

        self.out.fillBranch("JetID_Tight", jetid)
        self.out.fillBranch("nJetIDTight", sum(1 for x in jetid if x))
        return True


JetIDProducer_2024 = lambda: JetIDProducer()
JetIDProducer = lambda: JetIDProducer()

# Backward-compatible import name.  It no longer produces JetPUID_SF branches.
JPUIDSFProducer = lambda: JetIDProducer()
