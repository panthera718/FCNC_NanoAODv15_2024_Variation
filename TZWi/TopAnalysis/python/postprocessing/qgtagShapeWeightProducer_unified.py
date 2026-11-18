import os
import correctionlib
import ROOT
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class QGtagWeightProducer(Module):
    """
    Unified version of:
      - qgtagShapeWeightProducer.py
      - qgtagShapeWeightProducer_jetVarUp.py
      - qgtagShapeWeightProducer_jetVarDown.py
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
        self.json_path = os.path.join(cmssw_base, "src/TZWi/TopAnalysis/data/BTV/UL2018/qgtagging_v0.json")

        if not os.path.exists(self.json_path):
            raise RuntimeError("Q/G-tagging JSON file not found: %s" % self.json_path)

        print(f"Loading Q/G-tagging SF from correctionlib: {self.json_path}")

        self.cset = correctionlib.CorrectionSet.from_file(self.json_path)
        self.qg_corr = self.cset["deepJet_shape"]

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
        Dynamically maps keyword arguments to the exact positional arguments.
        Includes a fallback to 'central' (or 1.0) for missing systematics.
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
                args.append(kwargs.get("working_point", "M"))
            elif name == "flavor" or name == "absflavor":
                args.append(int(kwargs.get(name, 0)))
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
        if hasattr(self, "qg_corr"):
            del self.qg_corr
        if hasattr(self, "cset"):
            del self.cset

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch(self._branch_name("QGtagWeight"), "F")
        self.out.branch(self._branch_name("QGtagWeight_up"), "F")
        self.out.branch(self._branch_name("QGtagWeight_down"), "F")

    def analyze(self, event):
        tree = event._tree

        weight_central = 1.0
        weight_up = 1.0
        weight_down = 1.0

        n_good_jet = self._get_tree_value(tree, "nGoodJet", 0)

        if n_good_jet > 0:
            pts = self._get_tree_value(tree, "GoodJet_pt")
            etas = self._get_tree_value(tree, "GoodJet_eta")
            orig_idxs = self._get_tree_value(tree, "GoodJet_index")
            qg_scores = self._get_tree_value(tree, "GoodJet_DeepFlavQG")

            for jet_index in range(n_good_jet):
                jet_pt = pts[jet_index]
                jet_eta = etas[jet_index]
                abseta = abs(jet_eta)

                if jet_pt < 30.0 or abseta > 2.5:
                    continue

                orig_jet_idx = orig_idxs[jet_index]
                try:
                    parton_flavor = int(abs(event.Jet_partonFlavour[orig_jet_idx]))
                except (AttributeError, RuntimeError, IndexError, TypeError):
                    parton_flavor = int(abs(self._get_tree_value(tree, "GoodJet_hadronFlavour")[jet_index]))

                qvg = float(1.0 - float(qg_scores[jet_index]))

                eval_kwargs = {
                    "systematic": "central",
                    "pt": float(jet_pt),
                    "eta": float(abseta),
                    "abseta": float(abseta),
                    "flavor": parton_flavor,
                    "absflavor": parton_flavor,
                    "qvg": qvg,
                    "discr": qvg,
                }

                sf_central = self.evaluate_corr(self.qg_corr, **eval_kwargs)

                eval_kwargs["systematic"] = "up"
                sf_up = self.evaluate_corr(self.qg_corr, **eval_kwargs)

                eval_kwargs["systematic"] = "down"
                sf_down = self.evaluate_corr(self.qg_corr, **eval_kwargs)

                weight_central *= sf_central
                weight_up *= sf_up
                weight_down *= sf_down

        self.out.fillBranch(self._branch_name("QGtagWeight"), weight_central)
        self.out.fillBranch(self._branch_name("QGtagWeight_up"), weight_up)
        self.out.fillBranch(self._branch_name("QGtagWeight_down"), weight_down)

        return True


qgtagWeight = lambda: QGtagWeightProducer()


def _make_up_lambda(jet_var):
    return lambda: QGtagWeightProducer(jetVar=jet_var, jetSyst="Up")


def _make_down_lambda(jet_var):
    return lambda: QGtagWeightProducer(jetVar=jet_var, jetSyst="Down")


for _jet_var in QGtagWeightProducer.JET_VARS:
    globals()["qgtagWeightUp_%s" % _jet_var] = _make_up_lambda(_jet_var)
    globals()["qgtagWeightDown_%s" % _jet_var] = _make_down_lambda(_jet_var)
