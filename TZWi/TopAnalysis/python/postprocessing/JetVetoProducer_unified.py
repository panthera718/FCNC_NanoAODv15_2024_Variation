import os
import math
import correctionlib

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class JetVetoProducer(Module, object):
    """
    Run-3 2024 NanoAODv15 jet-veto-map producer using the official JME JSON.

    This replaces the old Run-2 ROOT histogram hot-jet map.  It must run after
    FCNCTriLepton, because the default evaluation uses the selected GoodJet
    collection.  The output convention keeps JetVetoSum for downstream
    compatibility:

      JetVetoSum = sum of non-zero veto-map values over selected GoodJets
      nJetVeto   = number of selected GoodJets in vetoed eta-phi regions

    The event-level veto should be applied as JetVetoSum == 0.
    """

    VALID_DIRECTIONS = {"central", "up", "down"}
    JET_VARS = [
        "AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias", "Fragmentation",
        "SinglePionECAL", "SinglePionHCAL", "FlavorQCD", "TimePtEta",
        "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB",
        "RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeBal",
        "RelativeSample", "RelativeFSR", "RelativeStatFSR", "RelativeStatEC",
        "RelativeStatHF", "PileUpDataMC", "PileUpPtRef", "PileUpPtBB",
        "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF", "FlavorPureQuark",
        "FlavorPureCharm", "FlavorPureBottom",
    ]

    def __init__(self, *args, **kwargs):
        self.jetVar = kwargs.get("jetVar")
        self.direction = kwargs.get("direction", "central")
        if self.direction not in self.VALID_DIRECTIONS:
            raise ValueError("direction must be one of %s" % sorted(self.VALID_DIRECTIONS))
        if self.direction != "central" and not self.jetVar:
            raise ValueError("Up/Down JetVetoProducer requires jetVar to be set")

        cmssw_base = os.environ.get("CMSSW_BASE", "")
        self.json_path = kwargs.get(
            "jsonPath",
            os.path.join(cmssw_base, "src/TZWi/TopAnalysis/data/JetVeto/2024/jetvetomaps.json.gz"),
        )
        self.correction_name = kwargs.get("correctionName", "Summer24Prompt24_RunBCDEFGHI_V1")
        self.map_type = kwargs.get("mapType", "jetvetomap")
        self.min_jet_pt = float(kwargs.get("minJetPt", 30.0))

        if not os.path.exists(self.json_path):
            raise RuntimeError("JME jet-veto-map JSON file not found: %s" % self.json_path)

        print("[JetVetoProducer] Loading official 2024 jet-veto map: %s" % self.json_path)
        self.cset = correctionlib.CorrectionSet.from_file(self.json_path)
        self.corr = self.cset[self.correction_name]

    def beginJob(self):
        pass

    def endJob(self):
        for attr in ("corr", "cset"):
            if hasattr(self, attr):
                delattr(self, attr)

    def _suffix(self):
        if self.direction == "central":
            return ""
        return "_%s%s" % (self.jetVar, "Up" if self.direction == "up" else "Down")

    def _branch_name(self, base):
        return base + self._suffix()

    @staticmethod
    def _get_tree_value(tree, name, default=None):
        return getattr(tree, "b_out_%s" % name, default)

    @staticmethod
    def _len(obj):
        try:
            return len(obj)
        except Exception:
            try:
                return obj.size()
            except Exception:
                return 0

    @staticmethod
    def _at(obj, idx, default=0.0):
        if obj is None:
            return default
        try:
            return obj[idx]
        except Exception:
            try:
                return obj.At(idx)
            except Exception:
                return default

    @staticmethod
    def _wrap_phi(phi):
        phi = float(phi)
        while phi >= math.pi:
            phi -= 2.0 * math.pi
        while phi < -math.pi:
            phi += 2.0 * math.pi
        return phi

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self._branch_name("JetVetoSum"), "F")
        self.out.branch(self._branch_name("JetVetoMax"), "F")
        self.out.branch(self._branch_name("nJetVeto"), "i")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def _goodjet_inputs_from_fcnc(self, event):
        tree = event._tree
        suffix = self._suffix()
        n_name = "nGoodJet" + suffix
        eta_name = "GoodJet_eta" + suffix
        phi_name = "GoodJet_phi" + suffix
        pt_name = "GoodJet_pt" + suffix

        n_good = int(self._get_tree_value(tree, n_name, 0) or 0)
        etas = self._get_tree_value(tree, eta_name, None)
        phis = self._get_tree_value(tree, phi_name, None)
        pts = self._get_tree_value(tree, pt_name, None)
        if n_good <= 0 or etas is None or phis is None:
            return []

        jets = []
        for idx in range(n_good):
            pt = float(self._at(pts, idx, self.min_jet_pt))
            if pt < self.min_jet_pt:
                continue
            jets.append((float(self._at(etas, idx, 99.0)), float(self._at(phis, idx, 99.0))))
        return jets

    def _raw_jet_inputs(self, event):
        jets = []
        try:
            for jet in Collection(event, "Jet"):
                if float(getattr(jet, "pt", 0.0)) < self.min_jet_pt:
                    continue
                jets.append((float(jet.eta), float(jet.phi)))
        except Exception:
            pass
        return jets

    def _veto_value(self, eta, phi):
        try:
            return float(self.corr.evaluate(self.map_type, float(eta), self._wrap_phi(phi)))
        except Exception:
            return 0.0

    def analyze(self, event):
        jets = self._goodjet_inputs_from_fcnc(event)
        if not jets:
            jets = self._raw_jet_inputs(event)

        veto_sum = 0.0
        veto_max = 0.0
        n_veto = 0
        for eta, phi in jets:
            value = self._veto_value(eta, phi)
            if value != 0.0:
                n_veto += 1
                veto_sum += value
                veto_max = max(veto_max, abs(value))

        tree = event._tree
        for base, value in (
            ("JetVetoSum", float(veto_sum)),
            ("JetVetoMax", float(veto_max)),
            ("nJetVeto", int(n_veto)),
        ):
            name = self._branch_name(base)
            setattr(tree, "b_out_%s" % name, value)
            self.out.fillBranch(name, value)

        return True


JVProducer = lambda: JetVetoProducer(direction="central")


def _make_lambda(direction, jetVar):
    return lambda direction=direction, jetVar=jetVar: JetVetoProducer(direction=direction, jetVar=jetVar)


for _jetVar in JetVetoProducer.JET_VARS:
    globals()["JVProducerUp_%s" % _jetVar] = _make_lambda("up", _jetVar)
    globals()["JVProducerDown_%s" % _jetVar] = _make_lambda("down", _jetVar)
