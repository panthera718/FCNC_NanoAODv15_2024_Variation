import os
import correctionlib
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


def is_relevant_syst_for_shape_corr(flavor_btv, syst):
    """
    Returns true if a flavor/syst combination is relevant.
    Note: correctionlib uses standard PDG ID flavors: 5(b), 4(c), 0(light)
    """
    if flavor_btv == 5:   # b-jet
        return syst in ["jes", "lf", "hfstats1", "hfstats2"]
    elif flavor_btv == 4: # c-jet
        return syst in ["cferr1", "cferr2"]
    elif flavor_btv == 0: # light-jet
        return syst in ["jes", "hf", "lfstats1", "lfstats2"]
    return False


class BtagShapeWeightProducer(Module):
    """
    Unified version of:
      - btagShapeWeightProducer.py
      - btagShapeWeightProducer_jetVarUp.py
      - btagShapeWeightProducer_jetVarDown.py

    Behavior is selected only by constructor arguments / lambda names.
    """

    JET_VARS = [
        "AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias",
        "Fragmentation", "SinglePionECAL", "SinglePionHCAL",
        "FlavorQCD",
        "TimePtEta", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
        "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF",
        "RelativeBal", "RelativeSample", "RelativeFSR", "RelativeStatFSR",
        "RelativeStatEC", "RelativeStatHF", "PileUpDataMC", "PileUpPtRef",
        "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF",
    ]

    def __init__(self, *args, **kwargs):
        self.jetVar = kwargs.get("jetVar")
        self.jetSyst = kwargs.get("jetSyst")
        if self.jetSyst not in (None, "Up", "Down"):
            self.jetSyst = None
        if self.jetSyst in ("Up", "Down") and not self.jetVar:
            self.jetSyst = None
            self.jetVar = None

        cmssw_base = os.environ.get("CMSSW_BASE", "")
        self.json_path = os.path.join(cmssw_base, "src/TZWi/TopAnalysis/data/BTV/UL2018/btagging.json.gz")

        if not os.path.exists(self.json_path):
            raise RuntimeError("BTV JSON file not found: %s" % self.json_path)

        print(f"Loading BTagCalibration from correctionlib: {self.json_path}")
        self.cset = correctionlib.CorrectionSet.from_file(self.json_path)
        self.btag_corr = self.cset["deepJet_shape"]

        self.systs = [
            "jes", "lf", "hf",
            "hfstats1", "hfstats2",
            "lfstats1", "lfstats2",
            "cferr1", "cferr2"
        ]

    def beginJob(self):
        pass

    def endJob(self):
        if hasattr(self, 'btag_corr'):
            del self.btag_corr
        if hasattr(self, 'cset'):
            del self.cset

    def _is_jetvar_mode(self):
        return self.jetVar is not None and self.jetSyst in ("Up", "Down")

    def _suffix(self):
        if not self._is_jetvar_mode():
            return ""
        return "_%s%s" % (self.jetVar, self.jetSyst)

    def _tree_attr(self, base):
        return "b_out_%s%s" % (base, self._suffix())

    def _branch_name(self, base):
        return "%s%s" % (base, self._suffix())

    def _get_tree_value(self, tree, base, default=None):
        return getattr(tree, self._tree_attr(base), default)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch(self._branch_name("BtagWeight"), "F")
        for syst in self.systs:
            self.out.branch(self._branch_name("BtagWeight_up_%s" % syst), "F")
            self.out.branch(self._branch_name("BtagWeight_down_%s" % syst), "F")

    def analyze(self, event):
        tree = event._tree
        weight_central = 1.0
        weights_up = {syst: 1.0 for syst in self.systs}
        weights_down = {syst: 1.0 for syst in self.systs}

        n_good_jet = self._get_tree_value(tree, "nGoodJet", 0)

        if n_good_jet > 0:
            pts = self._get_tree_value(tree, "GoodJet_pt")
            etas = self._get_tree_value(tree, "GoodJet_eta")
            flavs = self._get_tree_value(tree, "GoodJet_hadronFlavour")
            scores = self._get_tree_value(tree, "GoodJet_DeepFlavB")

            for jet_index in range(n_good_jet):
                jet_pt = pts[jet_index]
                jet_eta = etas[jet_index]

                if jet_pt < 25.0 or abs(jet_eta) > 2.4:
                    continue

                flav = abs(flavs[jet_index])
                flav_map = 0
                if flav == 5:
                    flav_map = 5
                elif flav == 4:
                    flav_map = 4

                score = float(scores[jet_index])

                sf_central = self.btag_corr.evaluate("central", flav_map, abs(jet_eta), jet_pt, score)
                if sf_central < 0.01:
                    sf_central = 1.0

                weight_central *= sf_central

                for syst in self.systs:
                    if is_relevant_syst_for_shape_corr(flav_map, syst):
                        sf_up = self.btag_corr.evaluate("up_" + syst, flav_map, abs(jet_eta), jet_pt, score)
                        sf_down = self.btag_corr.evaluate("down_" + syst, flav_map, abs(jet_eta), jet_pt, score)

                        if sf_up < 0.01:
                            sf_up = 1.0
                        if sf_down < 0.01:
                            sf_down = 1.0
                    else:
                        sf_up = sf_central
                        sf_down = sf_central

                    weights_up[syst] *= sf_up
                    weights_down[syst] *= sf_down

        self.out.fillBranch(self._branch_name("BtagWeight"), weight_central)
        for syst in self.systs:
            self.out.fillBranch(self._branch_name("BtagWeight_up_%s" % syst), weights_up[syst])
            self.out.fillBranch(self._branch_name("BtagWeight_down_%s" % syst), weights_down[syst])

        return True


btagWeight = lambda: BtagShapeWeightProducer()


def _make_up_lambda(jet_var):
    return lambda: BtagShapeWeightProducer(jetVar=jet_var, jetSyst="Up")


def _make_down_lambda(jet_var):
    return lambda: BtagShapeWeightProducer(jetVar=jet_var, jetSyst="Down")


for _jet_var in BtagShapeWeightProducer.JET_VARS:
    globals()["btagWeightUp_%s" % _jet_var] = _make_up_lambda(_jet_var)
    globals()["btagWeightDown_%s" % _jet_var] = _make_down_lambda(_jet_var)
