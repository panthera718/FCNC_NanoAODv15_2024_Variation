import math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class TriggerSFVariables2024(Module, object):
    """Write compact derived variables for 2D trigger-SF maps.

    This module must run after FCNCTriLepton and before the output branch
    selection is applied.  It reads the b_out_* cache filled by FCNCTriLepton
    and writes pT/|eta| variables that are stable across the four exact-3l
    channels.
    """

    BRANCHES = [
        "leadLep_pt", "subleadLep_pt", "trailLep_pt",
        "leadLep_absEta", "subleadLep_absEta", "trailLep_absEta",
        "leadEl_pt", "subleadEl_pt", "leadEl_absEta", "subleadEl_absEta",
        "leadMu_pt", "subleadMu_pt", "leadMu_absEta", "subleadMu_absEta",
    ]

    def __init__(self, *args, **kwargs):
        self.prefix = kwargs.get("prefix", "trigSF2D")
        self.out = None

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for name in self.BRANCHES:
            self.out.branch("%s_%s" % (self.prefix, name), "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    @staticmethod
    def _get_cached(tree, name, default=-999.0):
        return getattr(tree, "b_out_%s" % name, default)

    @staticmethod
    def _fill_float(out, tree, name, value):
        value = float(value)
        setattr(tree, "b_out_%s" % name, value)
        out.fillBranch(name, value)

    def analyze(self, event):
        tree = event._tree

        leptons = []
        electrons = []
        muons = []
        for idx in (1, 2, 3):
            pt = float(self._get_cached(tree, "Lepton%d_pt" % idx, -999.0))
            eta = float(self._get_cached(tree, "Lepton%d_eta" % idx, -999.0))
            pdg_id = int(self._get_cached(tree, "Lepton%d_pdgId" % idx, 0))
            item = {"pt": pt, "eta": eta, "absEta": abs(eta), "pdgId": pdg_id}
            if pt > -998.0:
                leptons.append(item)
                if abs(pdg_id) == 11:
                    electrons.append(item)
                elif abs(pdg_id) == 13:
                    muons.append(item)

        leptons.sort(key=lambda x: x["pt"], reverse=True)
        electrons.sort(key=lambda x: x["pt"], reverse=True)
        muons.sort(key=lambda x: x["pt"], reverse=True)

        values = {
            "leadLep_pt": -999.0,
            "subleadLep_pt": -999.0,
            "trailLep_pt": -999.0,
            "leadLep_absEta": -999.0,
            "subleadLep_absEta": -999.0,
            "trailLep_absEta": -999.0,
            "leadEl_pt": -999.0,
            "subleadEl_pt": -999.0,
            "leadEl_absEta": -999.0,
            "subleadEl_absEta": -999.0,
            "leadMu_pt": -999.0,
            "subleadMu_pt": -999.0,
            "leadMu_absEta": -999.0,
            "subleadMu_absEta": -999.0,
        }

        for rank, label in enumerate(["lead", "sublead", "trail"]):
            if len(leptons) > rank:
                values["%sLep_pt" % label] = leptons[rank]["pt"]
                values["%sLep_absEta" % label] = leptons[rank]["absEta"]

        if len(electrons) > 0:
            values["leadEl_pt"] = electrons[0]["pt"]
            values["leadEl_absEta"] = electrons[0]["absEta"]
        if len(electrons) > 1:
            values["subleadEl_pt"] = electrons[1]["pt"]
            values["subleadEl_absEta"] = electrons[1]["absEta"]
        if len(muons) > 0:
            values["leadMu_pt"] = muons[0]["pt"]
            values["leadMu_absEta"] = muons[0]["absEta"]
        if len(muons) > 1:
            values["subleadMu_pt"] = muons[1]["pt"]
            values["subleadMu_absEta"] = muons[1]["absEta"]

        for short_name in self.BRANCHES:
            branch_name = "%s_%s" % (self.prefix, short_name)
            self._fill_float(self.out, tree, branch_name, values[short_name])

        return True


trigSFVariables2024 = lambda: TriggerSFVariables2024()
