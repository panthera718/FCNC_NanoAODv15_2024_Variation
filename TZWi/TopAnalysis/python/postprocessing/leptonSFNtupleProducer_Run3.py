import math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class _ConstantArrayReader(object):
    def __init__(self, value=0.0, size=0):
        self.value = value
        self.size = int(size)

    def At(self, _):
        return self.value

    def GetSize(self):
        return self.size


class _ConstantValueReader(object):
    def __init__(self, value=0.0):
        self.value = value

    def Get(self):
        return self.value


class LeptonSFNtupleProducer(Module, object):
    """
    Run-3 NanoAODv15 promptMVA lepton-selection scale-factor NTuple producer.

    Event topology:
      - WZ-like exactly-three-relaxed-lepton event;
      - two strict, OS-SF tag leptons form the Z candidate;
      - the remaining relaxed lepton is the probe;
      - denominator = relaxed probe;
      - numerator   = the same probe also passes the strict promptMVA selection.

    Modes:
      - LepSFEl: electron probe
      - LepSFMu: muon probe

    The relaxed definitions are aligned with the analysis selection:
      - muon: mediumId, pfRelIso04_all < 0.40, |dxy| < 0.05,
              |dz| < 0.10, jetRelIso < 0.50;
      - electron: cutBased >= 3, convVeto, |dxy| < 0.05,
                  |dz| < 0.10, jetRelIso < 0.50, ECAL gap veto.

    Default 2024 strict promptMVA working points:
      - muon    : promptMVA > 0.70
      - electron: tightCharge >= 2 and promptMVA > 0.90
    """

    VALID_MODES = {"LepSFEl", "LepSFMu"}

    def __init__(self, *args, **kwargs):
        self.mode = kwargs.get("mode", "LepSFMu")
        self.year = str(kwargs.get("year", "2024"))
        if self.mode not in self.VALID_MODES:
            raise RuntimeError("Unsupported LeptonSFNtupleProducer mode: %s" % self.mode)

        self.z_mass = float(kwargs.get("zMass", 91.1876))
        self.z_window = float(kwargs.get("zWindow", 15.0))
        self.min_z_lead_pt = float(kwargs.get("minZLeadPt", 25.0))
        self.min_z_sublead_pt = float(kwargs.get("minZSubleadPt", 10.0))
        self.min_probe_pt = float(kwargs.get("minProbePt", 10.0))
        self.min_relaxed_pt = float(kwargs.get("minRelaxedPt", 10.0))
        self.min_trilepton_mass = float(kwargs.get("minTriLeptonMass", 100.0))
        self.require_exactly_three_relaxed = bool(kwargs.get("requireExactlyThreeRelaxed", True))
        self.require_trilepton_mass = bool(kwargs.get("requireTriLeptonMass", True))
        self.exclude_electron_ecal_gap = bool(kwargs.get("excludeElectronEcalGap", True))

        self.max_muon_abs_eta = float(kwargs.get("maxMuonAbsEta", 2.4))
        self.max_electron_abs_sc_eta = float(kwargs.get("maxElectronAbsScEta", 2.5))
        self.max_dxy = float(kwargs.get("maxDxy", 0.05))
        self.max_dz = float(kwargs.get("maxDz", 0.10))
        self.max_jet_rel_iso = float(kwargs.get("maxJetRelIso", 0.50))
        self.min_electron_cut_based = int(kwargs.get("minElectronCutBased", 3))
        self.max_muon_pf_rel_iso04 = float(kwargs.get("maxMuonPFRelIso04", 0.40))

        self.muon_prompt_mva_cut = float(kwargs.get("muonPromptMVACut", 0.70))
        self.electron_prompt_mva_cut = float(kwargs.get("electronPromptMVACut", 0.90))

        self.out = None

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        for attr in list(self.__dict__.keys()):
            if attr.startswith("b_"):
                delattr(self, attr)

    @staticmethod
    def _raw_tree(tree):
        return tree.tree() if hasattr(tree, "tree") and tree.tree() else tree

    def _has_branch(self, tree, branch_name):
        raw_tree = self._raw_tree(tree)
        return bool(hasattr(raw_tree, "GetBranch") and raw_tree.GetBranch(branch_name))

    def _array_reader_first(self, tree, branch_names, required=True, default=0.0):
        for branch_name in branch_names:
            if self._has_branch(tree, branch_name):
                return tree.arrayReader(branch_name)
        if required:
            raise RuntimeError("Missing required NanoAOD branch. Tried: %s" % ", ".join(branch_names))
        return _ConstantArrayReader(default)

    def _value_reader_first(self, tree, branch_names, required=True, default=0.0):
        for branch_name in branch_names:
            if self._has_branch(tree, branch_name):
                return tree.valueReader(branch_name)
        if required:
            raise RuntimeError("Missing required NanoAOD branch. Tried: %s" % ", ".join(branch_names))
        return _ConstantValueReader(default)

    @staticmethod
    def _uchar_to_int(value, default=0):
        """Convert NanoAOD UInt8/Bool-like PyROOT values to a Python int.

        In CMSSW/PyROOT, UInt8_t and Bool_t arrayReader values can be
        returned as one-character strings such as "\x01" instead of numeric
        Python integers.  Calling int("\x01") raises ValueError, and
        bool("\x00") is also True because a non-empty string is truthy.
        This helper must therefore be used for NanoAOD byte-sized branches
        such as *_genPartFlav, Electron_cutBased, Electron_tightCharge,
        Muon_mediumId, and Electron_convVeto.
        """
        if value is None:
            return int(default)
        if isinstance(value, bool):
            return int(value)
        if isinstance(value, bytes):
            if len(value) == 0:
                return int(default)
            return int(value[0])
        if isinstance(value, str):
            if len(value) == 0:
                return int(default)
            return ord(value[0])
        try:
            return int(value)
        except (TypeError, ValueError):
            try:
                return ord(value)
            except Exception:
                return int(default)

    @classmethod
    def _uchar_to_bool(cls, value, default=False):
        return bool(cls._uchar_to_int(value, int(default)))

    @staticmethod
    def _read_value(reader, default=0.0):
        if reader is None:
            return default
        try:
            value = reader.Get()
            try:
                return value[0]
            except TypeError:
                return value
        except Exception:
            return default

    @staticmethod
    def _p4(pt, eta, phi, mass):
        vec = ROOT.TLorentzVector()
        vec.SetPtEtaPhiM(float(pt), float(eta), float(phi), float(mass))
        return vec

    @staticmethod
    def _mt(lep_p4, met_pt, met_phi):
        dphi = ROOT.TVector2.Phi_mpi_pi(float(lep_p4.Phi()) - float(met_phi))
        return math.sqrt(max(0.0, 2.0 * float(lep_p4.Pt()) * float(met_pt) * (1.0 - math.cos(dphi))))

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        # Event identity and bookkeeping.  These are essential for deduplicating
        # data events across EGamma/Muon/MuonEG primary datasets in the SF step.
        self.out.branch("SF_run", "I")
        self.out.branch("SF_luminosityBlock", "I")
        self.out.branch("SF_event", "L")
        self.out.branch("SF_isData", "O")
        self.out.branch("SF_weight", "F")
        self.out.branch("SF_genWeight", "F")
        self.out.branch("SF_nRelaxedLepton", "I")
        self.out.branch("SF_nStrictLepton", "I")

        # Probe branches used by make_lepton_sf_promptMVA.py.
        self.out.branch("SF_probe_flavor", "I")       # 11 for electron, 13 for muon
        self.out.branch("SF_probe_pdgId", "I")
        self.out.branch("SF_probe_charge", "I")
        self.out.branch("SF_probe_pt", "F")
        self.out.branch("SF_probe_eta", "F")
        self.out.branch("SF_probe_absEta", "F")
        self.out.branch("SF_probe_scEta", "F")
        self.out.branch("SF_probe_absScEta", "F")
        self.out.branch("SF_probe_phi", "F")
        self.out.branch("SF_probe_mass", "F")
        self.out.branch("SF_probe_promptMVA", "F")
        self.out.branch("SF_probe_jetRelIso", "F")
        self.out.branch("SF_probe_pfRelIso04_all", "F")
        self.out.branch("SF_probe_cutBased", "I")
        self.out.branch("SF_probe_mediumId", "O")
        self.out.branch("SF_probe_tightCharge", "I")
        self.out.branch("SF_probe_convVeto", "O")
        self.out.branch("SF_probe_passRelax", "O")
        self.out.branch("SF_probe_passStrict", "O")
        self.out.branch("SF_probe_passPromptMVA", "O")
        self.out.branch("SF_probe_genPartFlav", "I")
        self.out.branch("SF_probe_isPromptMC", "O")

        # Tag/Z/WZ-event branches for validation plots and purity checks.
        self.out.branch("SF_tag1_pdgId", "I")
        self.out.branch("SF_tag2_pdgId", "I")
        self.out.branch("SF_tag1_pt", "F")
        self.out.branch("SF_tag2_pt", "F")
        self.out.branch("SF_tag1_eta", "F")
        self.out.branch("SF_tag2_eta", "F")
        self.out.branch("SF_tag1_phi", "F")
        self.out.branch("SF_tag2_phi", "F")
        self.out.branch("SF_tag1_promptMVA", "F")
        self.out.branch("SF_tag2_promptMVA", "F")
        self.out.branch("SF_tag1_jetRelIso", "F")
        self.out.branch("SF_tag2_jetRelIso", "F")
        self.out.branch("SF_tag1_pfRelIso04_all", "F")
        self.out.branch("SF_tag2_pfRelIso04_all", "F")
        self.out.branch("SF_tag1_cutBased", "I")
        self.out.branch("SF_tag2_cutBased", "I")
        self.out.branch("SF_tag1_mediumId", "O")
        self.out.branch("SF_tag2_mediumId", "O")
        self.out.branch("SF_tag1_tightCharge", "I")
        self.out.branch("SF_tag2_tightCharge", "I")
        self.out.branch("SF_z_mass", "F")
        self.out.branch("SF_z_pt", "F")
        self.out.branch("SF_z_eta", "F")
        self.out.branch("SF_z_phi", "F")
        self.out.branch("SF_trilepton_mass", "F")
        self.out.branch("SF_trilepton_pt", "F")
        self.out.branch("SF_w_mt", "F")
        self.out.branch("SF_met_pt", "F")
        self.out.branch("SF_met_phi", "F")

        # Scalar readers.
        self.b_run = self._value_reader_first(inputTree, ["run"], required=False, default=0)
        self.b_lumi = self._value_reader_first(inputTree, ["luminosityBlock"], required=False, default=0)
        self.b_event = self._value_reader_first(inputTree, ["event"], required=False, default=0)
        self.b_genWeight = self._value_reader_first(inputTree, ["genWeight"], required=False, default=1.0)
        self.has_genWeight = self._has_branch(inputTree, "genWeight")

        self.b_met_pt = self._value_reader_first(inputTree, ["PuppiMET_pt", "MET_T1_pt", "MET_pt"], required=False, default=0.0)
        self.b_met_phi = self._value_reader_first(inputTree, ["PuppiMET_phi", "MET_T1_phi", "MET_phi"], required=False, default=0.0)

        # Muon readers.
        self.b_nMuon = self._value_reader_first(inputTree, ["nMuon"], required=False, default=0)
        for var_name in ["pt", "eta", "phi", "mass", "charge", "dxy", "dz", "jetRelIso", "pfRelIso04_all", "promptMVA", "genPartFlav"]:
            setattr(
                self,
                "b_Muon_%s" % var_name,
                self._array_reader_first(inputTree, ["Muon_%s" % var_name], required=False, default=0.0),
            )
        self.b_Muon_mediumId = self._array_reader_first(inputTree, ["Muon_mediumId"], required=False, default=0)

        # Electron readers.
        self.b_nElectron = self._value_reader_first(inputTree, ["nElectron"], required=False, default=0)
        for var_name in ["pt", "eta", "deltaEtaSC", "phi", "mass", "charge", "dxy", "dz", "jetRelIso", "promptMVA", "genPartFlav"]:
            setattr(
                self,
                "b_Electron_%s" % var_name,
                self._array_reader_first(inputTree, ["Electron_%s" % var_name], required=False, default=0.0),
            )
        self.b_Electron_cutBased = self._array_reader_first(inputTree, ["Electron_cutBased"], required=False, default=0)
        self.b_Electron_tightCharge = self._array_reader_first(inputTree, ["Electron_tightCharge"], required=False, default=0)
        self.b_Electron_convVeto = self._array_reader_first(inputTree, ["Electron_convVeto"], required=False, default=1)

    def _build_muons(self):
        out = []
        n_muon = int(self._read_value(self.b_nMuon, 0))
        for i in range(n_muon):
            pt = float(self.b_Muon_pt.At(i))
            eta = float(self.b_Muon_eta.At(i))
            phi = float(self.b_Muon_phi.At(i))
            mass = float(self.b_Muon_mass.At(i))
            charge = int(self.b_Muon_charge.At(i))
            prompt_mva = float(self.b_Muon_promptMVA.At(i))
            jet_rel_iso = float(self.b_Muon_jetRelIso.At(i))
            pf_rel_iso04_all = float(self.b_Muon_pfRelIso04_all.At(i))
            medium_id = self._uchar_to_bool(self.b_Muon_mediumId.At(i))
            dxy = abs(float(self.b_Muon_dxy.At(i)))
            dz = abs(float(self.b_Muon_dz.At(i)))
            gen_part_flav = self._uchar_to_int(self.b_Muon_genPartFlav.At(i))

            pass_relax = (
                pt >= self.min_relaxed_pt and
                abs(eta) < self.max_muon_abs_eta and
                medium_id and
                pf_rel_iso04_all < self.max_muon_pf_rel_iso04 and
                dxy < self.max_dxy and
                dz < self.max_dz and
                jet_rel_iso < self.max_jet_rel_iso
            )
            pass_prompt_mva = prompt_mva > self.muon_prompt_mva_cut
            pass_strict = pass_relax and pass_prompt_mva

            p4 = self._p4(pt, eta, phi, mass)
            out.append({
                "kind": "mu",
                "flavor": 13,
                "index": i,
                "pt": pt,
                "eta": eta,
                "scEta": eta,
                "phi": phi,
                "mass": mass,
                "charge": charge,
                "pdgId": -charge * 13,
                "promptMVA": prompt_mva,
                "jetRelIso": jet_rel_iso,
                "pfRelIso04_all": pf_rel_iso04_all,
                "cutBased": -1,
                "mediumId": medium_id,
                "tightCharge": -1,
                "convVeto": True,
                "passRelax": pass_relax,
                "passPromptMVA": pass_prompt_mva,
                "passStrict": pass_strict,
                "genPartFlav": gen_part_flav,
                "isPromptMC": gen_part_flav in (1, 15),
                "p4": p4,
            })
        return out

    def _electron_sc_eta(self, i):
        return float(self.b_Electron_eta.At(i)) + float(self.b_Electron_deltaEtaSC.At(i))

    def _build_electrons(self):
        out = []
        n_electron = int(self._read_value(self.b_nElectron, 0))
        for i in range(n_electron):
            pt = float(self.b_Electron_pt.At(i))
            eta = float(self.b_Electron_eta.At(i))
            sc_eta = self._electron_sc_eta(i)
            abs_sc_eta = abs(sc_eta)
            phi = float(self.b_Electron_phi.At(i))
            mass = float(self.b_Electron_mass.At(i))
            charge = int(self.b_Electron_charge.At(i))
            prompt_mva = float(self.b_Electron_promptMVA.At(i))
            jet_rel_iso = float(self.b_Electron_jetRelIso.At(i))
            cut_based = self._uchar_to_int(self.b_Electron_cutBased.At(i))
            tight_charge = self._uchar_to_int(self.b_Electron_tightCharge.At(i))
            conv_veto = self._uchar_to_bool(self.b_Electron_convVeto.At(i), default=True)
            dxy = abs(float(self.b_Electron_dxy.At(i)))
            dz = abs(float(self.b_Electron_dz.At(i)))
            gen_part_flav = self._uchar_to_int(self.b_Electron_genPartFlav.At(i))

            in_acceptance = abs_sc_eta < self.max_electron_abs_sc_eta
            if self.exclude_electron_ecal_gap and 1.4442 < abs_sc_eta < 1.5660:
                in_acceptance = False

            pass_relax = (
                pt >= self.min_relaxed_pt and
                in_acceptance and
                cut_based >= self.min_electron_cut_based and
                conv_veto and
                dxy < self.max_dxy and
                dz < self.max_dz and
                jet_rel_iso < self.max_jet_rel_iso
            )
            pass_prompt_mva = prompt_mva > self.electron_prompt_mva_cut
            pass_strict = pass_relax and tight_charge >= 2 and pass_prompt_mva

            p4 = self._p4(pt, eta, phi, mass)
            out.append({
                "kind": "el",
                "flavor": 11,
                "index": i,
                "pt": pt,
                "eta": eta,
                "scEta": sc_eta,
                "phi": phi,
                "mass": mass,
                "charge": charge,
                "pdgId": -charge * 11,
                "promptMVA": prompt_mva,
                "jetRelIso": jet_rel_iso,
                "pfRelIso04_all": -1.0,
                "cutBased": cut_based,
                "mediumId": False,
                "tightCharge": tight_charge,
                "convVeto": conv_veto,
                "passRelax": pass_relax,
                "passPromptMVA": pass_prompt_mva,
                "passStrict": pass_strict,
                "genPartFlav": gen_part_flav,
                "isPromptMC": gen_part_flav in (1, 15),
                "p4": p4,
            })
        return out

    @staticmethod
    def _same_object(a, b):
        return a["kind"] == b["kind"] and a["index"] == b["index"]

    def _choose_candidate(self, relaxed):
        target_flavor = 11 if self.mode == "LepSFEl" else 13
        best = None
        strict = [lep for lep in relaxed if lep["passStrict"]]

        for i in range(len(strict)):
            for j in range(i + 1, len(strict)):
                tag1 = strict[i]
                tag2 = strict[j]
                if tag1["flavor"] != tag2["flavor"]:
                    continue
                if tag1["charge"] * tag2["charge"] >= 0:
                    continue

                tag_pts = sorted([tag1["pt"], tag2["pt"]], reverse=True)
                if tag_pts[0] < self.min_z_lead_pt or tag_pts[1] < self.min_z_sublead_pt:
                    continue

                z_p4 = tag1["p4"] + tag2["p4"]
                z_mass = float(z_p4.M())
                z_diff = abs(z_mass - self.z_mass)
                if z_diff >= self.z_window:
                    continue

                probes = [
                    lep for lep in relaxed
                    if lep["flavor"] == target_flavor and
                    lep["pt"] >= self.min_probe_pt and
                    (not self._same_object(lep, tag1)) and
                    (not self._same_object(lep, tag2))
                ]
                for probe in probes:
                    tri_p4 = tag1["p4"] + tag2["p4"] + probe["p4"]
                    tri_mass = float(tri_p4.M())
                    if self.require_trilepton_mass and tri_mass < self.min_trilepton_mass:
                        continue

                    # Prefer the candidate closest to the Z mass; if tied, prefer
                    # the highest-pT probe for determinism.
                    score = (z_diff, -probe["pt"])
                    if best is None or score < best["score"]:
                        best = {
                            "score": score,
                            "tag1": tag1,
                            "tag2": tag2,
                            "probe": probe,
                            "z_p4": z_p4,
                            "tri_p4": tri_p4,
                            "z_mass": z_mass,
                            "tri_mass": tri_mass,
                        }
        return best

    def analyze(self, event):
        muons = self._build_muons()
        electrons = self._build_electrons()
        relaxed = [lep for lep in (muons + electrons) if lep["passRelax"]]
        strict = [lep for lep in relaxed if lep["passStrict"]]

        if self.require_exactly_three_relaxed and len(relaxed) != 3:
            return False
        if len(relaxed) < 3:
            return False

        cand = self._choose_candidate(relaxed)
        if cand is None:
            return False

        tag1 = cand["tag1"]
        tag2 = cand["tag2"]
        probe = cand["probe"]
        z_p4 = cand["z_p4"]
        tri_p4 = cand["tri_p4"]

        met_pt = float(self._read_value(self.b_met_pt, 0.0))
        met_phi = float(self._read_value(self.b_met_phi, 0.0))
        w_mt = self._mt(probe["p4"], met_pt, met_phi)

        gen_weight = float(self._read_value(self.b_genWeight, 1.0))
        is_data = not self.has_genWeight
        weight = 1.0 if is_data else gen_weight

        self.out.fillBranch("SF_run", int(self._read_value(self.b_run, 0)))
        self.out.fillBranch("SF_luminosityBlock", int(self._read_value(self.b_lumi, 0)))
        self.out.fillBranch("SF_event", int(self._read_value(self.b_event, 0)))
        self.out.fillBranch("SF_isData", bool(is_data))
        self.out.fillBranch("SF_weight", float(weight))
        self.out.fillBranch("SF_genWeight", float(gen_weight))
        self.out.fillBranch("SF_nRelaxedLepton", int(len(relaxed)))
        self.out.fillBranch("SF_nStrictLepton", int(len(strict)))

        self.out.fillBranch("SF_probe_flavor", int(probe["flavor"]))
        self.out.fillBranch("SF_probe_pdgId", int(probe["pdgId"]))
        self.out.fillBranch("SF_probe_charge", int(probe["charge"]))
        self.out.fillBranch("SF_probe_pt", float(probe["pt"]))
        self.out.fillBranch("SF_probe_eta", float(probe["eta"]))
        self.out.fillBranch("SF_probe_absEta", abs(float(probe["eta"])))
        self.out.fillBranch("SF_probe_scEta", float(probe["scEta"]))
        self.out.fillBranch("SF_probe_absScEta", abs(float(probe["scEta"])))
        self.out.fillBranch("SF_probe_phi", float(probe["phi"]))
        self.out.fillBranch("SF_probe_mass", float(probe["mass"]))
        self.out.fillBranch("SF_probe_promptMVA", float(probe["promptMVA"]))
        self.out.fillBranch("SF_probe_jetRelIso", float(probe["jetRelIso"]))
        self.out.fillBranch("SF_probe_pfRelIso04_all", float(probe["pfRelIso04_all"]))
        self.out.fillBranch("SF_probe_cutBased", int(probe["cutBased"]))
        self.out.fillBranch("SF_probe_mediumId", bool(probe["mediumId"]))
        self.out.fillBranch("SF_probe_tightCharge", int(probe["tightCharge"]))
        self.out.fillBranch("SF_probe_convVeto", bool(probe["convVeto"]))
        self.out.fillBranch("SF_probe_passRelax", bool(probe["passRelax"]))
        self.out.fillBranch("SF_probe_passStrict", bool(probe["passStrict"]))
        self.out.fillBranch("SF_probe_passPromptMVA", bool(probe["passPromptMVA"]))
        self.out.fillBranch("SF_probe_genPartFlav", int(probe["genPartFlav"]))
        self.out.fillBranch("SF_probe_isPromptMC", bool(probe["isPromptMC"]))

        self.out.fillBranch("SF_tag1_pdgId", int(tag1["pdgId"]))
        self.out.fillBranch("SF_tag2_pdgId", int(tag2["pdgId"]))
        self.out.fillBranch("SF_tag1_pt", float(tag1["pt"]))
        self.out.fillBranch("SF_tag2_pt", float(tag2["pt"]))
        self.out.fillBranch("SF_tag1_eta", float(tag1["eta"]))
        self.out.fillBranch("SF_tag2_eta", float(tag2["eta"]))
        self.out.fillBranch("SF_tag1_phi", float(tag1["phi"]))
        self.out.fillBranch("SF_tag2_phi", float(tag2["phi"]))
        self.out.fillBranch("SF_tag1_promptMVA", float(tag1["promptMVA"]))
        self.out.fillBranch("SF_tag2_promptMVA", float(tag2["promptMVA"]))
        self.out.fillBranch("SF_tag1_jetRelIso", float(tag1["jetRelIso"]))
        self.out.fillBranch("SF_tag2_jetRelIso", float(tag2["jetRelIso"]))
        self.out.fillBranch("SF_tag1_pfRelIso04_all", float(tag1["pfRelIso04_all"]))
        self.out.fillBranch("SF_tag2_pfRelIso04_all", float(tag2["pfRelIso04_all"]))
        self.out.fillBranch("SF_tag1_cutBased", int(tag1["cutBased"]))
        self.out.fillBranch("SF_tag2_cutBased", int(tag2["cutBased"]))
        self.out.fillBranch("SF_tag1_mediumId", bool(tag1["mediumId"]))
        self.out.fillBranch("SF_tag2_mediumId", bool(tag2["mediumId"]))
        self.out.fillBranch("SF_tag1_tightCharge", int(tag1["tightCharge"]))
        self.out.fillBranch("SF_tag2_tightCharge", int(tag2["tightCharge"]))
        self.out.fillBranch("SF_z_mass", float(cand["z_mass"]))
        self.out.fillBranch("SF_z_pt", float(z_p4.Pt()))
        self.out.fillBranch("SF_z_eta", float(z_p4.Eta()))
        self.out.fillBranch("SF_z_phi", float(z_p4.Phi()))
        self.out.fillBranch("SF_trilepton_mass", float(cand["tri_mass"]))
        self.out.fillBranch("SF_trilepton_pt", float(tri_p4.Pt()))
        self.out.fillBranch("SF_w_mt", float(w_mt))
        self.out.fillBranch("SF_met_pt", float(met_pt))
        self.out.fillBranch("SF_met_phi", float(met_phi))

        return True


# Factories expected by nano_postproc.py -I module object.
leptonSF_LepSFEl_2024 = lambda: LeptonSFNtupleProducer(
    mode="LepSFEl",
    year="2024",
    minElectronCutBased=3,
    maxMuonPFRelIso04=0.40,
)
leptonSF_LepSFMu_2024 = lambda: LeptonSFNtupleProducer(
    mode="LepSFMu",
    year="2024",
    minElectronCutBased=3,
    maxMuonPFRelIso04=0.40,
)

# Aliases for scripts that build factory names without the leptonSF_ prefix.
LepSFEl_2024 = leptonSF_LepSFEl_2024
LepSFMu_2024 = leptonSF_LepSFMu_2024
