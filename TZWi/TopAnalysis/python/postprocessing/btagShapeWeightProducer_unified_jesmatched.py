import os
import correctionlib
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class BtagShapeWeightProducer(Module):
    """
    2024 UParTAK4 b-tagging scale-factor producer.

    The Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15 BTV JSON uses
    UParTAK4 working-point corrections:
      - UParTAK4_wp_values
      - UParTAK4_comb    for b/c jets
      - UParTAK4_mujets  for b/c jets in ttbar-dilepton-like regions
      - UParTAK4_light   for light jets

    This module keeps the historical output branch prefix BtagWeight for
    downstream compatibility, but the calculation is now WP-based and uses the
    GoodJet_DeepFlavB branch filled from Jet_btagUParTAK4B by the Run-3
    FCNCTriLepton producer.

    Important limitation:
      A mathematically complete WP event reweighting for nBjet categories also
      needs MC tagging-efficiency maps for untagged jets.  The official 2024
      btagging.json.gz supplied here does not contain those efficiencies.  The
      event weight produced here therefore applies the official WP SF to tagged
      selected jets.  This matches the information available in the JSON and is
      the safest drop-in replacement for the old shape-SF module without adding
      external efficiency inputs.
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
        "FlavorPureQuark", "FlavorPureCharm", "FlavorPureBottom",
    ]

    # Available in the attached 2024 UParTAK4 JSON for UParTAK4_comb.
    # UParTAK4_light only has central/up/down/correlated/uncorrelated; for the
    # detailed b/c-only sources light jets are kept at their central SF.
    SYSTS_2024 = [
        "up", "down",
        "up_correlated", "down_correlated",
        "up_uncorrelated", "down_uncorrelated",
        "up_jes", "down_jes",
        "up_jer", "down_jer",
        "up_pileup", "down_pileup",
        "up_statistic", "down_statistic",
        "up_mur", "down_mur",
        "up_muf", "down_muf",
        "up_pdfas", "down_pdfas",
        "up_isrdef", "down_isrdef",
        "up_fsrdef", "down_fsrdef",
        "up_topmass", "down_topmass",
        "up_bfragmentation", "down_bfragmentation",
        "up_hdamp", "down_hdamp",
        "up_type3", "down_type3",
    ]

    WP_FALLBACK_2024 = {
        "L": 0.0246,
        "M": 0.1272,
        "T": 0.4648,
        "XT": 0.6298,
        "XXT": 0.9739,
    }

    def __init__(self, *args, **kwargs):
        self.jetVar = kwargs.get("jetVar")
        self.jetSyst = kwargs.get("jetSyst")
        if self.jetSyst not in (None, "Up", "Down"):
            self.jetSyst = None
        if self.jetSyst in ("Up", "Down") and not self.jetVar:
            self.jetSyst = None
            self.jetVar = None

        self.working_point = str(kwargs.get("workingPoint", "L")).upper()
        if self.working_point not in self.WP_FALLBACK_2024:
            raise RuntimeError("Unsupported UParTAK4 working point: %s" % self.working_point)

        # Use 'comb' by default, as recommended for general regions.  Switch to
        # 'mujets' explicitly only for ttbar-dilepton-enriched analyses where the
        # comb calibration phase-space overlap is a concern.
        self.bc_calibration = str(kwargs.get("bcCalibration", "comb")).lower()
        if self.bc_calibration not in ("comb", "mujets"):
            raise RuntimeError("bcCalibration must be 'comb' or 'mujets', got: %s" % self.bc_calibration)

        self.min_btag_score = float(kwargs.get("btagWPValue", self.WP_FALLBACK_2024[self.working_point]))
        self.min_jet_pt = float(kwargs.get("minJetPt", 30.0))
        self.max_abs_eta = float(kwargs.get("maxAbsEta", 2.5))
        self.min_corr_pt = float(kwargs.get("minCorrectionPt", 20.0))
        self.max_corr_pt = float(kwargs.get("maxCorrectionPt", 1000.0))

        cmssw_base = os.environ.get("CMSSW_BASE", "")
        self.json_path = kwargs.get(
            "jsonPath",
            os.path.join(cmssw_base, "src/TZWi/TopAnalysis/data/BTV/2024/btagging.json.gz"),
        )

        if not os.path.exists(self.json_path):
            raise RuntimeError("BTV JSON file not found: %s" % self.json_path)

        print("Loading 2024 UParTAK4 b-tagging SF from correctionlib: %s" % self.json_path)
        self.cset = correctionlib.CorrectionSet.from_file(self.json_path)
        self.wp_corr = self.cset["UParTAK4_wp_values"]
        self.bc_corr = self.cset["UParTAK4_%s" % self.bc_calibration]
        self.light_corr = self.cset["UParTAK4_light"]

        # Prefer the JSON WP value over the hard-coded fallback.
        try:
            self.min_btag_score = float(self.wp_corr.evaluate(self.working_point))
        except Exception:
            pass

        self.bc_valid_systs = self._correction_systematics("UParTAK4_%s" % self.bc_calibration)
        self.light_valid_systs = self._correction_systematics("UParTAK4_light")

        # Output labels without the up_/down_ prefix duplication.  Example:
        #   up_jes -> BtagWeight_up_jes
        output_systs = []
        for label in self.SYSTS_2024:
            if label in ("up", "down"):
                continue
            if label.startswith("up_"):
                output_systs.append(label[3:])
            elif label.startswith("down_"):
                output_systs.append(label[5:])
        self.output_systs = sorted(set(output_systs))

    def beginJob(self):
        pass

    def endJob(self):
        for attr in ["light_corr", "bc_corr", "wp_corr", "cset"]:
            if hasattr(self, attr):
                delattr(self, attr)

    def _correction_systematics(self, correction_name):
        # correctionlib's Python object does not expose a stable public list of
        # category keys across all versions, so keep the authoritative list from
        # the JSON schema and validate by trial evaluation in _safe_sf.
        if correction_name == "UParTAK4_light":
            return {"central", "up", "down", "up_correlated", "down_correlated", "up_uncorrelated", "down_uncorrelated"}
        return {
            "central",
            "up", "down",
            "up_correlated", "down_correlated",
            "up_uncorrelated", "down_uncorrelated",
            "up_jes", "down_jes",
            "up_jer", "down_jer",
            "up_pileup", "down_pileup",
            "up_statistic", "down_statistic",
            "up_mur", "down_mur",
            "up_muf", "down_muf",
            "up_pdfas", "down_pdfas",
            "up_isrdef", "down_isrdef",
            "up_fsrdef", "down_fsrdef",
            "up_topmass", "down_topmass",
            "up_bfragmentation", "down_bfragmentation",
            "up_hdamp", "down_hdamp",
            "up_type3", "down_type3",
        }

    def _is_jetvar_mode(self):
        return self.jetVar is not None and self.jetSyst in ("Up", "Down")

    def _suffix(self):
        if not self._is_jetvar_mode():
            return ""
        return "_%s%s" % (self.jetVar, self.jetSyst)

    def _branch_name(self, base):
        return "%s%s" % (base, self._suffix())

    def _active_jes_branch_base(self):
        if not self._is_jetvar_mode():
            return None
        if self.jetSyst == "Up":
            return "BtagWeight_up_jes"
        return "BtagWeight_down_jes"

    def _active_jes_label(self):
        if not self._is_jetvar_mode():
            return None
        return "up_jes" if self.jetSyst == "Up" else "down_jes"

    @staticmethod
    def _get_tree_value(tree, name, default=None):
        return getattr(tree, "b_out_%s" % name, default)

    @staticmethod
    def _len(obj):
        try:
            return len(obj)
        except Exception:
            if hasattr(obj, "size"):
                return obj.size()
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

    def _flavor_from_goodjet_or_raw(self, event, goodjet_index, local_index):
        tree = event._tree
        good_flavs = self._get_tree_value(tree, "GoodJet_hadronFlavour", None)
        if good_flavs is not None:
            return int(abs(self._at(good_flavs, local_index, 0)))

        try:
            jets = Collection(event, "Jet")
            jet = jets[int(goodjet_index)]
            if not hasattr(jet, "hadronFlavour"):
                return None
            return int(abs(getattr(jet, "hadronFlavour")))
        except Exception:
            return None

    @staticmethod
    def _flavor_for_btv(flavor):
        flavor = abs(int(flavor))
        if flavor == 5:
            return 5
        if flavor == 4:
            return 4
        return 0

    def _clamp_for_correction(self, pt, abs_eta):
        pt = float(pt)
        abs_eta = abs(float(abs_eta))

        # The JSON has flow handling in pt, but eta flow is an error.  Clamp only
        # for numerical safety; jets outside the valid eta acceptance are skipped
        # before calling this function.
        eps_pt = 1.0e-4
        eps_eta = 1.0e-6
        if pt < self.min_corr_pt:
            pt = self.min_corr_pt + eps_pt
        if pt >= self.max_corr_pt:
            pt = self.max_corr_pt - eps_pt
        if abs_eta < 0.0:
            abs_eta = 0.0
        if abs_eta >= self.max_abs_eta:
            abs_eta = self.max_abs_eta - eps_eta
        return pt, abs_eta

    def _safe_sf(self, systematic, flavor_btv, pt, abs_eta):
        if abs_eta >= self.max_abs_eta:
            return 1.0
        if pt < self.min_jet_pt:
            return 1.0

        corr = self.bc_corr if flavor_btv in (4, 5) else self.light_corr
        valid = self.bc_valid_systs if flavor_btv in (4, 5) else self.light_valid_systs
        label = systematic if systematic in valid else "central"
        pt_eval, eta_eval = self._clamp_for_correction(pt, abs_eta)

        try:
            sf = corr.evaluate(label, self.working_point, int(flavor_btv), eta_eval, pt_eval)
        except Exception:
            if label != "central":
                try:
                    sf = corr.evaluate("central", self.working_point, int(flavor_btv), eta_eval, pt_eval)
                except Exception:
                    sf = 1.0
            else:
                sf = 1.0

        sf = float(sf)
        if not (0.01 <= sf <= 10.0):
            return 1.0
        return sf

    def _event_weight(self, event, systematic):
        tree = event._tree
        n_good_jet = int(self._get_tree_value(tree, "nGoodJet", 0) or 0)
        pts = self._get_tree_value(tree, "GoodJet_pt", None)
        etas = self._get_tree_value(tree, "GoodJet_eta", None)
        scores = self._get_tree_value(tree, "GoodJet_DeepFlavB", None)
        indices = self._get_tree_value(tree, "GoodJet_index", None)

        if n_good_jet <= 0 or pts is None or etas is None or scores is None:
            return 1.0

        weight = 1.0
        for local_index in range(n_good_jet):
            pt = float(self._at(pts, local_index, 0.0))
            eta = float(self._at(etas, local_index, 99.0))
            score = float(self._at(scores, local_index, -999.0))
            if pt < self.min_jet_pt or abs(eta) >= self.max_abs_eta:
                continue
            if score <= self.min_btag_score:
                continue

            goodjet_index = int(self._at(indices, local_index, local_index))
            raw_flavor = self._flavor_from_goodjet_or_raw(event, goodjet_index, local_index)
            if raw_flavor is None:
                # Data or flavor-less private samples: b-tag SF is MC-only.
                continue
            flavor_btv = self._flavor_for_btv(raw_flavor)
            weight *= self._safe_sf(systematic, flavor_btv, pt, abs(eta))

        return float(weight)

    @staticmethod
    def _label_for_output(direction, syst):
        return "%s_%s" % (direction, syst)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        if self._is_jetvar_mode():
            self.out.branch(self._branch_name(self._active_jes_branch_base()), "F")
        else:
            self.out.branch(self._branch_name("BtagWeight"), "F")
            self.out.branch(self._branch_name("BtagWeight_up"), "F")
            self.out.branch(self._branch_name("BtagWeight_down"), "F")
            for syst in self.output_systs:
                self.out.branch(self._branch_name("BtagWeight_up_%s" % syst), "F")
                self.out.branch(self._branch_name("BtagWeight_down_%s" % syst), "F")

    def analyze(self, event):
        if self._is_jetvar_mode():
            self.out.fillBranch(
                self._branch_name(self._active_jes_branch_base()),
                self._event_weight(event, self._active_jes_label()),
            )
            return True

        self.out.fillBranch(self._branch_name("BtagWeight"), self._event_weight(event, "central"))
        self.out.fillBranch(self._branch_name("BtagWeight_up"), self._event_weight(event, "up"))
        self.out.fillBranch(self._branch_name("BtagWeight_down"), self._event_weight(event, "down"))
        for syst in self.output_systs:
            self.out.fillBranch(
                self._branch_name("BtagWeight_up_%s" % syst),
                self._event_weight(event, self._label_for_output("up", syst)),
            )
            self.out.fillBranch(
                self._branch_name("BtagWeight_down_%s" % syst),
                self._event_weight(event, self._label_for_output("down", syst)),
            )
        return True


btagWeight = lambda: BtagShapeWeightProducer()


def _make_up_lambda(jet_var):
    return lambda: BtagShapeWeightProducer(jetVar=jet_var, jetSyst="Up")


def _make_down_lambda(jet_var):
    return lambda: BtagShapeWeightProducer(jetVar=jet_var, jetSyst="Down")


for _jet_var in BtagShapeWeightProducer.JET_VARS:
    globals()["btagWeightUp_%s" % _jet_var] = _make_up_lambda(_jet_var)
    globals()["btagWeightDown_%s" % _jet_var] = _make_down_lambda(_jet_var)
