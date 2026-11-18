import os
import correctionlib
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class CtagShapeWeightProducer(Module):
    """
    Unified version of:
      - ctagShapeWeightProducer.py
      - ctagShapeWeightProducer_jetVarUp.py
      - ctagShapeWeightProducer_jetVarDown.py
    """

    JET_VARS = [
        "AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias",
        "Fragmentation", "SinglePionECAL", "SinglePionHCAL",
        "FlavorQCD",
        "TimePtEta", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
        "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF",
        "RelativeBal", "RelativeSample", "RelativeFSR",
        "RelativeStatFSR", "RelativeStatEC", "RelativeStatHF",
        "PileUpDataMC", "PileUpPtRef", "PileUpPtBB",
        "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF",
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
        self.json_path = os.path.join(cmssw_base, "src/TZWi/TopAnalysis/data/BTV/UL2018/ctagging.json.gz")

        if not os.path.exists(self.json_path):
            raise RuntimeError("C-tagging JSON file not found: %s" % self.json_path)

        self.systs = [
            "Extrap", "Interp", "Stat",
            "XSec_BRUnc_DYJets_b", "XSec_BRUnc_DYJets_c", "XSec_BRUnc_WJets_c",
            "jer", "jesTotal", "puWeight"
        ]

        print(f"Loading C-tagging shape SF from correctionlib: {self.json_path}")
        self.cset = correctionlib.CorrectionSet.from_file(self.json_path)
        self.ctag_corr = self.cset["deepJet_shape"]

    def _is_jetvar_mode(self):
        return self.jetVar is not None and self.jetSyst in ("Up", "Down")

    def _suffix(self):
        if not self._is_jetvar_mode():
            return ""
        return "_%s%s" % (self.jetVar, self.jetSyst)

    def _branch_name(self, base):
        return "%s%s" % (base, self._suffix())

    def _tree_attr(self, base):
        return "b_out_%s%s" % (base, self._suffix())

    def _get_tree_value(self, tree, base, default=None):
        if self._is_jetvar_mode():
            return getattr(tree, self._tree_attr(base), default)
        return getattr(tree, "b_out_%s" % base, default)

    def evaluate_corr(self, corr, **kwargs):
        """
        Dynamically map keyword arguments to the exact positional arguments
        required by the specific correctionlib JSON schema.
        Fallback to 'central' (or 1.0) for missing systematics.
        """
        args = []
        syst_index = -1

        for i, inp in enumerate(corr.inputs):
            name = inp.name.lower()
            if name == "systematic":
                syst_index = i

            if name in kwargs:
                args.append(kwargs[name])
            elif name == "working_point":
                args.append(kwargs.get("working_point", "shape"))
            else:
                args.append(kwargs.get(inp.name, 0.0))

        try:
            return corr.evaluate(*args)
        except Exception:
            if syst_index >= 0 and args[syst_index] != "central":
                args[syst_index] = "central"
                try:
                    return corr.evaluate(*args)
                except Exception:
                    return 1.0
            return 1.0

    def beginJob(self):
        pass

    def endJob(self):
        if hasattr(self, "ctag_corr"):
            del self.ctag_corr
        if hasattr(self, "cset"):
            del self.cset

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch(self._branch_name("CtagWeight"), "F")
        for syst in self.systs:
            self.out.branch(self._branch_name("CtagWeight_up_%s" % syst), "F")
            self.out.branch(self._branch_name("CtagWeight_down_%s" % syst), "F")

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
            cvls = self._get_tree_value(tree, "GoodJet_DeepFlavCvL")
            cvbs = self._get_tree_value(tree, "GoodJet_DeepFlavCvB")

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

                cvl = float(cvls[jet_index])
                cvb = float(cvbs[jet_index])

                eval_kwargs = {
                    "systematic": "central",
                    "flavor": flav_map,
                    "cvl": cvl,
                    "cvb": cvb,
                    "pt": jet_pt,
                    "eta": abs(jet_eta),
                    "abseta": abs(jet_eta),
                }

                sf_central = self.evaluate_corr(self.ctag_corr, **eval_kwargs)
                weight_central *= sf_central

                for syst in self.systs:
                    eval_kwargs["systematic"] = "up_" + syst
                    sf_up = self.evaluate_corr(self.ctag_corr, **eval_kwargs)

                    eval_kwargs["systematic"] = "down_" + syst
                    sf_down = self.evaluate_corr(self.ctag_corr, **eval_kwargs)

                    weights_up[syst] *= sf_up
                    weights_down[syst] *= sf_down

        self.out.fillBranch(self._branch_name("CtagWeight"), weight_central)
        for syst in self.systs:
            self.out.fillBranch(self._branch_name("CtagWeight_up_%s" % syst), weights_up[syst])
            self.out.fillBranch(self._branch_name("CtagWeight_down_%s" % syst), weights_down[syst])

        return True


ctagWeight = lambda: CtagShapeWeightProducer()


def _make_up_lambda(jet_var):
    return lambda: CtagShapeWeightProducer(jetVar=jet_var, jetSyst="Up")


def _make_down_lambda(jet_var):
    return lambda: CtagShapeWeightProducer(jetVar=jet_var, jetSyst="Down")


for _jet_var in CtagShapeWeightProducer.JET_VARS:
    globals()["ctagWeightUp_%s" % _jet_var] = _make_up_lambda(_jet_var)
    globals()["ctagWeightDown_%s" % _jet_var] = _make_down_lambda(_jet_var)
