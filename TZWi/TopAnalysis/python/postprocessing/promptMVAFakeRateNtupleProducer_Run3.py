import math
import os
import re

try:
    import correctionlib
except Exception:
    correctionlib = None

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


def _env_float(name, default):
    try:
        return float(os.environ.get(name, default))
    except Exception:
        return float(default)


def _env_int(name, default):
    try:
        return int(os.environ.get(name, default))
    except Exception:
        return int(default)


def _env_optional_float(name):
    value = os.environ.get(name, "")
    if value is None or str(value).strip() == "":
        return None
    try:
        return float(value)
    except Exception:
        return None


class PromptMVAFakeRateNtupleProducer(Module, object):
    """
    Run-3 NanoAODv15 promptMVA fake-rate NTuple producer for 2024.

    Nominal measurement strategy:
      * input data stream: Run2024.JetMET
      * trigger: probe-lepton-independent PFJet trigger OR from 2024_FR_JetMET_Trigger.yaml
      * FREl: one loose/fakeable electron probe in a QCD-enriched one-lepton region
      * FRMu: one loose/fakeable muon probe in a QCD-enriched one-lepton region

    Denominator:
      loose/fakeable probe lepton.

    Numerator:
      the same probe lepton passing the strict promptMVA selection.

    Stored regions:
      FRRegion=1: denominator event with loose-not-tight probe.
      FRRegion=2: denominator event with tight probe, also used as numerator.
    """

    FR_DEN = 1
    FR_NUM = 2

    # Nominal JetMET fake-rate trigger paths. The YAML file can override this list.
    DEFAULT_JETMET_TRIGGERS = (
        "HLT_PFJet200",
        "HLT_PFJet260",
        "HLT_PFJet320",
        "HLT_PFJet400",
        "HLT_PFJet450",
        "HLT_PFJet500",
        "HLT_PFJet550",
    )

    # 2024 Summer24 UParTAK4 working points from the BTV btagging JSON.
    # The main trilepton analysis uses the loose (L) WP for nBjet.
    UPARTAK4_BTAG_WP_2024 = {
        "L": 0.0246,
        "M": 0.1272,
        "T": 0.4648,
        "XT": 0.6298,
        "XXT": 0.9739,
    }
    BTAG_WP_CODE = {"L": 0, "M": 1, "T": 2, "XT": 3, "XXT": 4}

    def __init__(self, *args, **kwargs):
        self.mode = kwargs.get("mode", "FRMu")
        self.year = str(kwargs.get("year", "2024"))
        self.require_trigger = bool(kwargs.get("require_trigger", True))
        self.apply_met_filters = bool(kwargs.get("apply_met_filters", True))

        self.min_probe_pt = float(kwargs.get("min_probe_pt", _env_float("FR_MIN_PROBE_PT", 10.0)))
        self.min_lead_jet_pt = float(kwargs.get("min_lead_jet_pt", _env_float("FR_MIN_LEAD_JET_PT", 250.0)))
        self.max_met_plus_mt = float(kwargs.get("max_met_plus_mt", _env_float("FR_MAX_MET_PLUS_MT", 50.0)))
        self.max_loose_leptons = int(kwargs.get("max_loose_leptons", _env_int("FR_MAX_LOOSE_LEPTONS", 1)))

        # Analysis-dependent loose/tight lepton choices.
        self.mu_promptmva_wp = float(kwargs.get("mu_promptmva_wp", _env_float("FR_MU_PROMPTMVA_WP", 0.70)))
        self.el_promptmva_wp = float(kwargs.get("el_promptmva_wp", _env_float("FR_EL_PROMPTMVA_WP", 0.90)))
        self.max_jet_rel_iso = float(kwargs.get("max_jet_rel_iso", _env_float("FR_MAX_JET_REL_ISO", 0.50)))
        self.max_mu_mini_iso = float(kwargs.get("max_mu_mini_iso", _env_float("FR_MAX_MU_MINI_ISO", 0.40)))
        self.min_el_cutbased = int(kwargs.get("min_el_cutbased", _env_int("FR_MIN_EL_CUTBASED", 3)))

        # Keep the PFJet-trigger jet and the analysis/b-tag jets conceptually
        # separate.  The high-pT jet only establishes the trigger plateau.
        self.min_analysis_jet_pt = float(
            kwargs.get("min_analysis_jet_pt", _env_float("FR_ANALYSIS_MIN_JET_PT", 30.0))
        )
        self.max_analysis_jet_abs_eta = float(
            kwargs.get("max_analysis_jet_abs_eta", _env_float("FR_ANALYSIS_MAX_JET_ABS_ETA", 4.7))
        )
        self.min_btag_jet_pt = float(
            kwargs.get("min_btag_jet_pt", _env_float("FR_BTAG_MIN_JET_PT", 20.0))
        )
        self.max_btag_jet_abs_eta = float(
            kwargs.get("max_btag_jet_abs_eta", _env_float("FR_BTAG_MAX_JET_ABS_ETA", 2.5))
        )
        self.jet_lepton_dr = float(
            kwargs.get("jet_lepton_dr", _env_float("FR_JET_LEPTON_DR", 0.4))
        )

        self.btag_wp = str(kwargs.get("btag_wp", os.environ.get("FR_BTAG_WP", "L"))).upper()
        if self.btag_wp not in self.UPARTAK4_BTAG_WP_2024:
            raise RuntimeError(
                "Unsupported 2024 UParTAK4 working point %s; choose one of %s"
                % (self.btag_wp, ", ".join(self.UPARTAK4_BTAG_WP_2024))
            )
        env_wp_value = _env_optional_float("FR_BTAG_WP_VALUE")
        default_wp_value = self.UPARTAK4_BTAG_WP_2024[self.btag_wp]
        self.btag_wp_value = float(
            kwargs.get("btag_wp_value", default_wp_value if env_wp_value is None else env_wp_value)
        )

        # Evaluate the same official 2024 AK4PUPPI Tight jet ID as the main
        # trilepton producer.  A composition-based fallback is retained.
        self.apply_official_jet_id = bool(
            kwargs.get("apply_official_jet_id", _env_int("FR_APPLY_OFFICIAL_JET_ID", 1))
        )
        self.jet_id_working_point = str(
            kwargs.get("jet_id_working_point", os.environ.get("FR_JETID_WP", "AK4PUPPI_Tight"))
        )
        cmssw_base = os.environ.get("CMSSW_BASE", "")
        self.jet_id_json_path = kwargs.get(
            "jet_id_json_path",
            os.environ.get(
                "FR_JETID_JSON",
                os.path.join(cmssw_base, "src/TZWi/TopAnalysis/data/JetID/2024/jetid.json.gz"),
            ),
        )
        self.jet_id_cset = None
        self.jet_id_corr = None
        if self.apply_official_jet_id:
            if correctionlib is None:
                print("[PromptMVAFakeRateNtupleProducer] correctionlib unavailable; using jet-ID fallback.")
            elif not os.path.exists(self.jet_id_json_path):
                print("[PromptMVAFakeRateNtupleProducer] jet-ID JSON not found: %s" % self.jet_id_json_path)
                print("[PromptMVAFakeRateNtupleProducer] using jet-ID fallback.")
            else:
                try:
                    self.jet_id_cset = correctionlib.CorrectionSet.from_file(self.jet_id_json_path)
                    self.jet_id_corr = self.jet_id_cset[self.jet_id_working_point]
                except Exception as exc:
                    self.jet_id_cset = None
                    self.jet_id_corr = None
                    print("[PromptMVAFakeRateNtupleProducer] failed to load jet-ID JSON: %s" % exc)
                    print("[PromptMVAFakeRateNtupleProducer] using jet-ID fallback.")

        self.trigger_yaml = kwargs.get("trigger_yaml", os.environ.get("FR_TRIGGER_YAML", ""))
        self.trigger_key = kwargs.get("trigger_key", os.environ.get("FR_TRIGGER_KEY", "Run2024.JetMET_FR"))
        self.jetmet_triggers = tuple(kwargs.get("jetmet_triggers", self.DEFAULT_JETMET_TRIGGERS))
        yaml_paths = self._paths_from_trigger_yaml(self.trigger_yaml, self.trigger_key)
        if yaml_paths:
            self.jetmet_triggers = tuple(yaml_paths)

        self.out = None

        if self.mode not in ("FRMu", "FREl"):
            raise RuntimeError("Unsupported PromptMVAFakeRateNtupleProducer mode: %s" % self.mode)


    def beginJob(self):
        print("[PromptMVAFakeRateNtupleProducer] mode=%s year=%s" % (self.mode, self.year))
        print("[PromptMVAFakeRateNtupleProducer] require_trigger=%s trigger_key=%s trigger_yaml=%s" % (
            self.require_trigger, self.trigger_key, self.trigger_yaml if self.trigger_yaml else "<not set>"))
        print("[PromptMVAFakeRateNtupleProducer] trigger paths: %s" % ", ".join(self.jetmet_triggers))
        print("[PromptMVAFakeRateNtupleProducer] min_probe_pt=%.1f min_lead_jet_pt=%.1f max_met_plus_mt=%.1f" % (
            self.min_probe_pt, self.min_lead_jet_pt, self.max_met_plus_mt))
        print(
            "[PromptMVAFakeRateNtupleProducer] analysis jets: pT>%.1f |eta|<%.1f; "
            "b jets: pT>%.1f |eta|<%.1f UParTAK4 %s WP=%.4f"
            % (
                self.min_analysis_jet_pt,
                self.max_analysis_jet_abs_eta,
                self.min_btag_jet_pt,
                self.max_btag_jet_abs_eta,
                self.btag_wp,
                self.btag_wp_value,
            )
        )
        print(
            "[PromptMVAFakeRateNtupleProducer] jet ID=%s JSON=%s loaded=%s"
            % (
                self.jet_id_working_point,
                self.jet_id_json_path,
                bool(self.jet_id_corr is not None),
            )
        )

    def endJob(self):
        self.jet_id_corr = None
        self.jet_id_cset = None

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for name, typ in [
            ("FRRegion", "I"),
            ("FRTrigger_pass", "O"),
            ("nLooseLepton", "I"),
            ("nTightLepton", "I"),
            ("nFailLep", "I"),
            ("frWeightSign", "I"),
            ("eventWeight", "F"),
            ("run", "I"),
            ("luminosityBlock", "I"),
            ("event", "L"),
            ("PuppiMET_pt_out", "F"),
            ("PuppiMET_phi_out", "F"),
            ("FR_nAnalysisJet", "I"),
            ("FR_nBJet", "I"),
            ("FR_btagCategory0GE1", "I"),
            ("FR_btagCategory01GE2", "I"),
            ("FR_btagWPCode", "I"),
            ("FR_btagWPValue", "F"),
            ("FR_maxBtagUParTAK4B", "F"),
            ("FR_leadBJetPt", "F"),
            ("FRProbe_assocJetUParTAK4B", "F"),
            ("FRProbe_assocJetIsBTagged", "O"),
            ("FRProbe_flavor", "I"),
            ("FRProbe_pt", "F"),
            ("FRProbe_eta", "F"),
            ("FRProbe_absEta", "F"),
            ("FRProbe_phi", "F"),
            ("FRProbe_mass", "F"),
            ("FRProbe_charge", "I"),
            ("FRProbe_promptMVA", "F"),
            ("FRProbe_miniPFRelIso_all", "F"),
            ("FRProbe_jetRelIso", "F"),
            ("FRProbe_jetPtRatio", "F"),
            ("FRProbe_jetPtRelv2", "F"),
            ("FRProbe_jetDF", "F"),
            ("FRProbe_jetIdx", "I"),
            ("FRProbe_isLoose", "O"),
            ("FRProbe_isTight", "O"),
            ("FRProbe_mt", "F"),
            ("FRProbe_met_plus_mt", "F"),
            ("FRProbe_leadJetPt", "F"),
            ("FRProbe_genPartFlav", "I"),
        ]:
            self.out.branch(name, typ)

        raw_tree = inputTree.tree() if hasattr(inputTree, "tree") and inputTree.tree() else inputTree
        if not (hasattr(raw_tree, "GetBranch") and raw_tree.GetBranch("Jet_btagUParTAK4B")):
            raise RuntimeError(
                "Missing required NanoAODv15 branch Jet_btagUParTAK4B; "
                "cannot build UParTAK4 b-jet categories."
            )

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    @staticmethod
    def _get(obj, name, default=0):
        """Safely read an attribute from a NanoAODTools Object/Event.

        NanoAODTools can raise RuntimeError, not AttributeError, when the
        corresponding branch is absent from the input tree. This is expected
        for MC-only branches such as Electron_genPartFlav and Muon_genPartFlav
        when running on data. In that case, return the requested default.
        """
        if not isinstance(name, str):
            return default
        try:
            return getattr(obj, name)
        except (AttributeError, RuntimeError, TypeError):
            return default

    @staticmethod
    def _delta_phi(phi1, phi2):
        dphi = phi1 - phi2
        while dphi > math.pi:
            dphi -= 2.0 * math.pi
        while dphi <= -math.pi:
            dphi += 2.0 * math.pi
        return dphi

    @classmethod
    def _delta_r(cls, eta1, phi1, eta2, phi2):
        return math.hypot(eta1 - eta2, cls._delta_phi(phi1, phi2))

    @classmethod
    def _mt(cls, pt, phi, met_pt, met_phi):
        return math.sqrt(max(0.0, 2.0 * pt * met_pt * (1.0 - math.cos(cls._delta_phi(phi, met_phi)))))

    @staticmethod
    def _read_trigger_value_from_yaml(path, key):
        if not path or not os.path.exists(path):
            return None
        try:
            import yaml
            with open(path) as f:
                data = yaml.safe_load(f) or {}
            value = data.get(key, None)
            if value is not None:
                return str(value)
        except Exception:
            pass

        # Fallback parser for simple YAML blocks like:
        # Run2024.JetMET_FR: " ... "
        try:
            text = open(path).read()
            pattern = r"(?ms)^\s*" + re.escape(key) + r"\s*:\s*\"(.*?)\""
            match = re.search(pattern, text)
            if match:
                return match.group(1)
        except Exception:
            pass
        return None

    @classmethod
    def _paths_from_trigger_yaml(cls, path, key):
        expr = cls._read_trigger_value_from_yaml(path, key)
        if not expr:
            return []
        # Extract HLT branch names. For the nominal FR key the expression is a pure OR.
        # This deliberately does not try to evaluate complex lepton-stream veto logic;
        # fake-rate production should use the simple JetMET_FR key.
        paths = re.findall(r"HLT_[A-Za-z0-9_]+", expr)
        seen = set()
        out = []
        for p in paths:
            if p not in seen:
                seen.add(p)
                out.append(p)
        return out

    def _passes_any_trigger(self, event, paths):
        if not self.require_trigger:
            return True

        found_any_configured_branch = False
        for path in paths:
            value = self._get(event, path, None)
            if value is None:
                continue
            found_any_configured_branch = True
            if bool(value):
                return True

        # Fail closed. If none of the configured trigger branches exists, this
        # catches YAML/branch-name mistakes immediately. If they exist but are
        # all false, the event simply fails the trigger requirement.
        return False

    def _passes_met_filters(self, event):
        if not self.apply_met_filters:
            return True
        filters = (
            "Flag_goodVertices",
            "Flag_globalSuperTightHalo2016Filter",
            "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_BadPFMuonFilter",
            "Flag_BadPFMuonDzFilter",
            "Flag_hfNoisyHitsFilter",
            "Flag_eeBadScFilter",
            "Flag_ecalBadCalibFilter",
        )
        for flag in filters:
            value = self._get(event, flag, None)
            # Some flags can be absent in particular NanoAOD campaigns. Missing
            # flags are not used to reject the event; present false flags reject it.
            if value is not None and not bool(value):
                return False
        return True

    def _jet_pt_ratio(self, lep):
        jet_rel_iso = float(self._get(lep, "jetRelIso", -1.0))
        if jet_rel_iso < -0.99:
            return -1.0
        return 1.0 / (1.0 + jet_rel_iso)

    def _is_muon_loose(self, mu):
        return (
            float(mu.pt) > self.min_probe_pt and
            abs(float(mu.eta)) < 2.4 and
            bool(self._get(mu, "mediumId", False)) and
            abs(float(self._get(mu, "dxy", 999.0))) < 0.05 and
            abs(float(self._get(mu, "dz", 999.0))) < 0.10 and
            float(self._get(mu, "miniPFRelIso_all", 999.0)) < self.max_mu_mini_iso and
            float(self._get(mu, "jetRelIso", 999.0)) < self.max_jet_rel_iso
        )

    def _is_muon_tight(self, mu):
        return self._is_muon_loose(mu) and float(self._get(mu, "promptMVA", -999.0)) > self.mu_promptmva_wp

    def _is_electron_loose(self, el):
        sc_eta = float(el.eta) + float(self._get(el, "deltaEtaSC", 0.0))
        in_gap = 1.4442 < abs(sc_eta) < 1.566
        return (
            float(el.pt) > self.min_probe_pt and
            abs(float(el.eta)) < 2.5 and
            not in_gap and
            int(self._get(el, "cutBased", 0)) >= self.min_el_cutbased and
            abs(float(self._get(el, "dxy", 999.0))) < 0.05 and
            abs(float(self._get(el, "dz", 999.0))) < 0.10 and
            float(self._get(el, "jetRelIso", 999.0)) < self.max_jet_rel_iso and
            bool(self._get(el, "convVeto", False))
        )

    def _is_electron_tight(self, el):
        return (
            self._is_electron_loose(el) and
            float(self._get(el, "promptMVA", -999.0)) > self.el_promptmva_wp and
            int(self._get(el, "tightCharge", 0)) >= 2
        )

    def _make_lepton_state(self, flavor, obj, idx):
        if flavor == 13:
            is_loose = self._is_muon_loose(obj)
            is_tight = self._is_muon_tight(obj)
        else:
            is_loose = self._is_electron_loose(obj)
            is_tight = self._is_electron_tight(obj)
        return {
            "flavor": int(flavor),
            "pt": float(obj.pt),
            "eta": float(obj.eta),
            "absEta": abs(float(obj.eta)),
            "phi": float(obj.phi),
            "mass": float(self._get(obj, "mass", 0.0)),
            "charge": int(self._get(obj, "charge", 0)),
            "promptMVA": float(self._get(obj, "promptMVA", -999.0)),
            "miniPFRelIso_all": float(self._get(obj, "miniPFRelIso_all", -1.0)),
            "jetRelIso": float(self._get(obj, "jetRelIso", -1.0)),
            "jetPtRatio": float(self._jet_pt_ratio(obj)),
            "jetPtRelv2": float(self._get(obj, "jetPtRelv2", -1.0)),
            "jetDF": float(self._get(obj, "jetDF", 0.0)),
            "jetIdx": int(self._get(obj, "jetIdx", -1)),
            "genPartFlav": int(self._get(obj, "genPartFlav", -1)),
            "isLoose": bool(is_loose),
            "isTight": bool(is_tight),
        }

    def _all_leptons(self, event):
        leptons = []
        for i, mu in enumerate(Collection(event, "Muon")):
            leptons.append(self._make_lepton_state(13, mu, i))
        for i, el in enumerate(Collection(event, "Electron")):
            leptons.append(self._make_lepton_state(11, el, i))
        leptons.sort(key=lambda x: x["pt"], reverse=True)
        return leptons

    @staticmethod
    def _to_int(value):
        if isinstance(value, str):
            return ord(value[0]) if value else 0
        if isinstance(value, bytes):
            return value[0] if value else 0
        try:
            return int(value)
        except Exception:
            return 0

    def _passes_official_ak4puppi_tight_jet_id(self, jet):
        if self.jet_id_corr is None:
            return None
        try:
            ch_mult = self._to_int(self._get(jet, "chMultiplicity", 0))
            ne_mult = self._to_int(self._get(jet, "neMultiplicity", 0))
            multiplicity = ch_mult + ne_mult
            value = self.jet_id_corr.evaluate(
                float(jet.eta),
                float(self._get(jet, "chHEF", 0.0)),
                float(self._get(jet, "neHEF", 0.0)),
                float(self._get(jet, "chEmEF", 0.0)),
                float(self._get(jet, "neEmEF", 0.0)),
                float(self._get(jet, "muEF", 0.0)),
                int(ch_mult),
                int(ne_mult),
                int(multiplicity),
            )
            return float(value) > 0.5
        except Exception:
            return None

    def _passes_tight_jet_id_from_composition(self, jet):
        # Same fallback used by the 2024 trilepton analysis producer.
        eta = abs(float(jet.eta))
        ne_hef = float(self._get(jet, "neHEF", 0.0))
        ne_emef = float(self._get(jet, "neEmEF", 0.0))
        ch_hef = float(self._get(jet, "chHEF", 0.0))
        ch_emef = float(self._get(jet, "chEmEF", 0.0))
        mu_ef = float(self._get(jet, "muEF", 0.0))
        ch_mult = self._to_int(self._get(jet, "chMultiplicity", 0))
        ne_mult = self._to_int(self._get(jet, "neMultiplicity", 0))
        n_const = self._to_int(self._get(jet, "nConstituents", ch_mult + ne_mult))

        if eta <= 2.6:
            return (
                ne_hef < 0.99 and
                ne_emef < 0.90 and
                n_const > 1 and
                mu_ef < 0.80 and
                ch_hef > 0.01 and
                ch_mult > 0 and
                ch_emef < 0.80
            )
        if eta <= 2.7:
            return ne_hef < 0.90 and ne_emef < 0.99 and mu_ef < 0.80
        if eta <= 3.0:
            return ne_emef < 0.99 and ne_mult > 2
        return ne_hef < 0.90 and ne_emef < 0.90 and ne_mult > 10

    def _passes_tight_jet_id(self, jet):
        official = self._passes_official_ak4puppi_tight_jet_id(jet)
        if official is not None:
            return bool(official)

        jet_id = self._get(jet, "jetId", None)
        if jet_id is not None:
            return int(jet_id) >= 2

        return self._passes_tight_jet_id_from_composition(jet)

    def _jet_is_clean(self, jet, leptons):
        jet_eta = float(jet.eta)
        jet_phi = float(jet.phi)
        for lep in leptons:
            if self._delta_r(jet_eta, jet_phi, lep["eta"], lep["phi"]) < self.jet_lepton_dr:
                return False
        return True

    def _trigger_jets(self, event, probe):
        # High-pT jets used only for the JetMET trigger plateau requirement.
        jets = []
        for jet in Collection(event, "Jet"):
            if float(jet.pt) <= self.min_lead_jet_pt:
                continue
            if abs(float(jet.eta)) >= self.max_analysis_jet_abs_eta:
                continue
            if not self._passes_tight_jet_id(jet):
                continue
            if not self._jet_is_clean(jet, [probe]):
                continue
            jets.append(jet)
        jets.sort(key=lambda x: float(x.pt), reverse=True)
        return jets

    def _count_analysis_jets_and_bjets(self, event, leptons):
        # Match the main 2024 trilepton analysis:
        #   analysis jets: pT > 30 GeV, |eta| < 4.7
        #   b-tag count : pT > 20 GeV, |eta| < 2.5
        # Both use AK4PUPPI Tight jet ID and DeltaR(jet, selected lepton) >= 0.4.
        n_analysis_jet = 0
        n_bjet = 0
        max_btag_score = -999.0
        lead_bjet_pt = -999.0

        for jet in Collection(event, "Jet"):
            if not self._passes_tight_jet_id(jet):
                continue
            if not self._jet_is_clean(jet, leptons):
                continue

            pt = float(jet.pt)
            abs_eta = abs(float(jet.eta))

            if pt > self.min_analysis_jet_pt and abs_eta < self.max_analysis_jet_abs_eta:
                n_analysis_jet += 1

            if not (pt > self.min_btag_jet_pt and abs_eta < self.max_btag_jet_abs_eta):
                continue

            score = self._get(jet, "btagUParTAK4B", None)
            if score is None:
                raise RuntimeError(
                    "Jet_btagUParTAK4B disappeared while processing the event; "
                    "the b-tag category cannot be defined."
                )
            score = float(score)
            max_btag_score = max(max_btag_score, score)
            if score > self.btag_wp_value:
                n_bjet += 1
                lead_bjet_pt = max(lead_bjet_pt, pt)

        return {
            "nAnalysisJet": int(n_analysis_jet),
            "nBJet": int(n_bjet),
            "category0GE1": 0 if n_bjet == 0 else 1,
            "category01GE2": min(int(n_bjet), 2),
            "maxBtagScore": float(max_btag_score),
            "leadBJetPt": float(lead_bjet_pt),
        }

    def _probe_associated_jet_btag(self, event, probe):
        # Diagnostic only.  FR_nBJet follows the analysis lepton-cleaned jet
        # definition, while this value records the probe-associated source jet.
        jet_idx = int(probe.get("jetIdx", -1))
        if jet_idx < 0:
            return -999.0, False
        jets = Collection(event, "Jet")
        if jet_idx >= len(jets):
            return -999.0, False
        score = self._get(jets[jet_idx], "btagUParTAK4B", None)
        if score is None:
            return -999.0, False
        score = float(score)
        return score, bool(score > self.btag_wp_value)

    def _fill_event(self, event, region, trigger_pass, probe, n_loose, n_tight, mt, met_pt, met_phi, lead_jet_pt, jet_summary, assoc_jet_btag):
        self.out.fillBranch("FRRegion", int(region))
        self.out.fillBranch("FRTrigger_pass", bool(trigger_pass))
        self.out.fillBranch("nLooseLepton", int(n_loose))
        self.out.fillBranch("nTightLepton", int(n_tight))
        self.out.fillBranch("nFailLep", int(n_loose - n_tight))
        self.out.fillBranch("frWeightSign", 1)
        gen_weight = float(self._get(event, "genWeight", 1.0))
        self.out.fillBranch("eventWeight", 1.0 if gen_weight >= 0.0 else -1.0)
        self.out.fillBranch("run", int(event.run))
        self.out.fillBranch("luminosityBlock", int(event.luminosityBlock))
        self.out.fillBranch("event", int(event.event))
        self.out.fillBranch("PuppiMET_pt_out", float(met_pt))
        self.out.fillBranch("PuppiMET_phi_out", float(met_phi))

        self.out.fillBranch("FR_nAnalysisJet", int(jet_summary["nAnalysisJet"]))
        self.out.fillBranch("FR_nBJet", int(jet_summary["nBJet"]))
        self.out.fillBranch("FR_btagCategory0GE1", int(jet_summary["category0GE1"]))
        self.out.fillBranch("FR_btagCategory01GE2", int(jet_summary["category01GE2"]))
        self.out.fillBranch("FR_btagWPCode", int(self.BTAG_WP_CODE[self.btag_wp]))
        self.out.fillBranch("FR_btagWPValue", float(self.btag_wp_value))
        self.out.fillBranch("FR_maxBtagUParTAK4B", float(jet_summary["maxBtagScore"]))
        self.out.fillBranch("FR_leadBJetPt", float(jet_summary["leadBJetPt"]))
        self.out.fillBranch("FRProbe_assocJetUParTAK4B", float(assoc_jet_btag[0]))
        self.out.fillBranch("FRProbe_assocJetIsBTagged", bool(assoc_jet_btag[1]))

        self.out.fillBranch("FRProbe_flavor", int(probe["flavor"]))
        self.out.fillBranch("FRProbe_pt", float(probe["pt"]))
        self.out.fillBranch("FRProbe_eta", float(probe["eta"]))
        self.out.fillBranch("FRProbe_absEta", float(probe["absEta"]))
        self.out.fillBranch("FRProbe_phi", float(probe["phi"]))
        self.out.fillBranch("FRProbe_mass", float(probe["mass"]))
        self.out.fillBranch("FRProbe_charge", int(probe["charge"]))
        self.out.fillBranch("FRProbe_promptMVA", float(probe["promptMVA"]))
        self.out.fillBranch("FRProbe_miniPFRelIso_all", float(probe["miniPFRelIso_all"]))
        self.out.fillBranch("FRProbe_jetRelIso", float(probe["jetRelIso"]))
        self.out.fillBranch("FRProbe_jetPtRatio", float(probe["jetPtRatio"]))
        self.out.fillBranch("FRProbe_jetPtRelv2", float(probe["jetPtRelv2"]))
        self.out.fillBranch("FRProbe_jetDF", float(probe["jetDF"]))
        self.out.fillBranch("FRProbe_jetIdx", int(probe["jetIdx"]))
        self.out.fillBranch("FRProbe_isLoose", bool(probe["isLoose"]))
        self.out.fillBranch("FRProbe_isTight", bool(probe["isTight"]))
        self.out.fillBranch("FRProbe_mt", float(mt))
        self.out.fillBranch("FRProbe_met_plus_mt", float(met_pt + mt))
        self.out.fillBranch("FRProbe_leadJetPt", float(lead_jet_pt))
        self.out.fillBranch("FRProbe_genPartFlav", int(probe["genPartFlav"]))

    def analyze(self, event):
        if not self._passes_met_filters(event):
            return False

        trigger_pass = self._passes_any_trigger(event, self.jetmet_triggers)
        if not trigger_pass:
            return False

        target_flavor = 13 if self.mode == "FRMu" else 11
        leptons = self._all_leptons(event)
        loose_leptons = [lep for lep in leptons if lep["isLoose"]]
        target_loose = [lep for lep in loose_leptons if lep["flavor"] == target_flavor]
        if len(target_loose) != 1:
            return False
        if len(loose_leptons) > self.max_loose_leptons:
            return False

        probe = target_loose[0]
        met_pt_value = self._get(event, "PuppiMET_pt", None)
        met_phi_value = self._get(event, "PuppiMET_phi", None)
        if met_pt_value is None:
            met_pt_value = self._get(event, "MET_pt", 0.0)
        if met_phi_value is None:
            met_phi_value = self._get(event, "MET_phi", 0.0)
        met_pt = float(met_pt_value)
        met_phi = float(met_phi_value)
        mt = self._mt(probe["pt"], probe["phi"], met_pt, met_phi)
        if met_pt + mt >= self.max_met_plus_mt:
            return False

        trigger_jets = self._trigger_jets(event, probe)
        if len(trigger_jets) < 1:
            return False
        lead_jet_pt = float(trigger_jets[0].pt)

        # nBJet is counted independently of the pT>250 GeV trigger jet.
        # Use the analysis-consistent pT/eta acceptance and UParTAK4 loose WP.
        jet_summary = self._count_analysis_jets_and_bjets(event, loose_leptons)
        assoc_jet_btag = self._probe_associated_jet_btag(event, probe)

        n_loose = len(loose_leptons)
        n_tight = sum(1 for lep in loose_leptons if lep["isTight"])
        region = self.FR_NUM if probe["isTight"] else self.FR_DEN
        self._fill_event(
            event,
            region,
            trigger_pass,
            probe,
            n_loose,
            n_tight,
            mt,
            met_pt,
            met_phi,
            lead_jet_pt,
            jet_summary,
            assoc_jet_btag,
        )
        return True


FRProducer_FRMu_2024 = lambda: PromptMVAFakeRateNtupleProducer(mode="FRMu", year="2024", require_trigger=True)
FRProducer_FREl_2024 = lambda: PromptMVAFakeRateNtupleProducer(mode="FREl", year="2024", require_trigger=True)

# Debug-only factories. Do not use for nominal fake-rate production.
FRProducerNoHLT_FRMu_2024 = lambda: PromptMVAFakeRateNtupleProducer(mode="FRMu", year="2024", require_trigger=False)
FRProducerNoHLT_FREl_2024 = lambda: PromptMVAFakeRateNtupleProducer(mode="FREl", year="2024", require_trigger=False)
