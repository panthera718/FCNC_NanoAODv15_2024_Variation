import math
import itertools
import os

try:
    import correctionlib
except Exception:
    correctionlib = None

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection  # kept for interface compatibility
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class _ConstantArrayReader(object):
    def __init__(self, value=0.0, size=0):
        self.value = value
        self.size = int(size)

    def At(self, _):
        return self.value

    def GetSize(self):
        return self.size


class _JetPtRatioFromRelIsoReader(object):
    def __init__(self, rel_iso_reader):
        self.rel_iso_reader = rel_iso_reader

    def At(self, i):
        rel_iso = float(self.rel_iso_reader.At(i))
        if rel_iso < 0.0:
            return 1.0
        return 1.0 / (1.0 + rel_iso)


class _JetComponentFromIndexReader(object):
    def __init__(self, lep_jet_idx_reader, jet_component_reader, default=0.0):
        self.lep_jet_idx_reader = lep_jet_idx_reader
        self.jet_component_reader = jet_component_reader
        self.default = default

    def At(self, i):
        idx = int(self.lep_jet_idx_reader.At(i))
        if idx < 0 or idx >= self.jet_component_reader.GetSize():
            return self.default
        return self.jet_component_reader.At(idx)


class _OptionalBranchMixin(object):
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

    def _value_reader_first(self, tree, branch_names, required=True):
        for branch_name in branch_names:
            if self._has_branch(tree, branch_name):
                return tree.valueReader(branch_name)
        if required:
            raise RuntimeError("Missing required NanoAOD branch. Tried: %s" % ", ".join(branch_names))
        return None

    @staticmethod
    def _event_value_first(event, branch_names):
        """Return the first available scalar event branch without using hasattr(event, ...).

        NanoAODTools' Event.__getattr__ calls readBranch(name).  Therefore
        hasattr(event, missing_branch) raises RuntimeError instead of returning
        False for missing NanoAOD branches, e.g. MET_T1_pt in Run-3
        NanoAODv15.  Check the underlying TTree branch list first, then read
        only branches that actually exist.
        """
        tree = getattr(event, "_tree", None)
        raw_tree = None
        if tree is not None:
            try:
                raw_tree = tree.tree() if hasattr(tree, "tree") and tree.tree() else tree
            except Exception:
                raw_tree = tree

        for branch_name in branch_names:
            has_branch = False
            try:
                has_branch = bool(
                    raw_tree is not None and
                    hasattr(raw_tree, "GetBranch") and
                    raw_tree.GetBranch(branch_name)
                )
            except Exception:
                has_branch = False

            if not has_branch:
                continue

            if tree is not None and hasattr(tree, "readBranch"):
                return tree.readBranch(branch_name)
            return getattr(event, branch_name)

        raise RuntimeError("Missing required event branch. Tried: %s" % ", ".join(branch_names))



class FCNCTriLepton(_OptionalBranchMixin, Module, object):
    """
    Pure-Python replacement for the previous
      - fcncTriLepton.py / central_fcncTriLepton.py
      - FCNCTriLeptonCppWorker.h / ForBTag_FCNCTriLeptonCppWorker.h
      - FCNCTriLeptonCppWorker.cc / ForBTag_FCNCTriLeptonCppWorker.cc

    The event-selection and branch-filling logic is ported from the original
    Python+C++ implementation so that this module can run standalone without
    loading the C++ worker.
    """

    VALID_MODES = {
        "MuElEl",
        "ElMuMu",
        "ElElEl",
        "MuMuMu",
        "TaTaTa",
        "ElTaTa",
        "MuTaTa",
        "TaElEl",
        "TaMuMu",
        "NPLMuElEl",
        "NPLElMuMu",
        "NPLElElEl",
        "NPLMuMuMu",
        "NPLTaTaTa",
        "NPLElTaTa",
        "NPLMuTaTa",
        "NPLTaElEl",
        "NPLTaMuMu",
    }

    # 2024 Summer24 UParTAK4 b-tag working points from BTV btagging.json.gz.
    # L/M/T/XT/XXT correspond to 10/1/0.1/0.05/0.01% light-jet mis-id rates.
    UParTAK4_BTAG_WP_2024 = {
        "L": 0.0246,
        "M": 0.1272,
        "T": 0.4648,
        "XT": 0.6298,
        "XXT": 0.9739,
    }

    def __init__(self, *args, **kwargs):
        self.mode = kwargs.get("mode")
        self.jetVar = kwargs.get("jetVar")
        self.jetSyst = kwargs.get("jetSyst")
        if self.jetSyst not in (None, "Up", "Down"):
            print('Jet systematic direction "%s" is not available.' % self.jetSyst)
            self.jetSyst = None
        if self.jetSyst in ("Up", "Down") and not self.jetVar:
            print('jetSyst="%s" requires jetVar to be set.' % self.jetSyst)
            self.jetSyst = None
            self.jetVar = None

        self.use_hadronFlavour = bool(
            kwargs.get("use_hadronFlavour", False) or
            (self.jetVar is not None and self.jetSyst in ("Up", "Down"))
        )
        if self.mode not in self.VALID_MODES:
            print('Mode name "%s" is not available.' % self.mode)
            self.mode = None

        self.minCjetCvLDiscr_ = 0.099
        self.minCjetCvBDiscr_ = 0.325
        self.minQjetQGDiscr_ = 0.27

        # 2024 uses the UParTAK4 AK4 b tagger.  Keep the downstream branch names
        # GoodJet_DeepFlav* for compatibility, but use Jet_btagUParTAK4* as the
        # preferred input branch and the 2024 BTV WP values for the tag decision.
        self.btagWP_ = str(kwargs.get("btagWP", "L")).upper()
        default_btag_wp = self.UParTAK4_BTAG_WP_2024.get(self.btagWP_, self.UParTAK4_BTAG_WP_2024["L"])
        self.minBjetBDiscr_ = float(kwargs.get("btagLooseWP", default_btag_wp))
        self.maxBTagJetAbsEta_ = float(kwargs.get("maxBTagJetAbsEta", 2.5))

        # Jet selection.  AN2025/031 uses Run-3 PUPPI AK4 jets and does not apply
        # the old CHS pileup-jet-ID requirement.  NanoAODv15 2024 does not store
        # Jet_puId or Jet_jetId for these AK4PUPPI jets; jet ID is evaluated from
        # the official JME jetid.json.gz, with a composition-based fallback only.
        self.minJetPt_ = float(kwargs.get("minJetPt", 30.0))
        self.maxJetAbsEta_ = float(kwargs.get("maxJetAbsEta", 4.7))
        self.applyPileupJetID_ = False

        # Official Run-3 AK4PUPPI tight jet ID from JME jetid.json.gz.
        # If the JSON or correctionlib is unavailable, fall back to the local
        # composition-based implementation using NanoAODv15 jet variables.
        self.applyOfficialJetID_ = bool(kwargs.get("applyOfficialJetID", True))
        self.jetIDWorkingPoint_ = str(kwargs.get("jetIDWorkingPoint", "AK4PUPPI_Tight"))
        cmssw_base = os.environ.get("CMSSW_BASE", "")
        self.jetIDJsonPath_ = kwargs.get(
            "jetIDJsonPath",
            os.path.join(cmssw_base, "src/TZWi/TopAnalysis/data/JetID/2024/jetid.json.gz"),
        )
        self.jetIDCorr_ = None
        if self.applyOfficialJetID_:
            if correctionlib is None:
                print("[FCNCTriLepton] correctionlib is not available; using composition fallback for jet ID.")
            elif not os.path.exists(self.jetIDJsonPath_):
                print("[FCNCTriLepton] JME jet-ID JSON not found: %s" % self.jetIDJsonPath_)
                print("[FCNCTriLepton] Using composition fallback for jet ID.")
            else:
                try:
                    self.jetIDCSet_ = correctionlib.CorrectionSet.from_file(self.jetIDJsonPath_)
                    self.jetIDCorr_ = self.jetIDCSet_[self.jetIDWorkingPoint_]
                    print("[FCNCTriLepton] Loaded official 2024 jet ID %s from %s" % (self.jetIDWorkingPoint_, self.jetIDJsonPath_))
                except Exception as exc:
                    self.jetIDCorr_ = None
                    print("[FCNCTriLepton] Failed to load official jet ID JSON: %s" % str(exc))
                    print("[FCNCTriLepton] Using composition fallback for jet ID.")

        # Run-3 prompt-lepton MVA working points.
        # AN2025/031 Table 8: strict muon promptMVA > 0.7.
        # AN2025/031 Table 9: strict electron promptMVA > 0.9 for 2024 NanoAODv15
        #                  and > 0.3 for 2022-2023 NanoAODv12.
        self.minMuonPromptMVA_ = float(kwargs.get("muonPromptMVACut", 0.7))
        self.minElectronPromptMVA_ = float(kwargs.get("electronPromptMVACut", 0.9))

        self._ttreereaderversion = -1
        self.out = None
    def beginJob(self):
        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        # No C++ worker remains, but keep the cleanup pattern for tree-attached
        # temporary outputs so the module stays friendly to the old workflow.
        for attr in list(self.__dict__.keys()):
            if attr.startswith("b_"):
                delattr(self, attr)

        if inputTree and hasattr(inputTree, 'tree') and inputTree.tree():
            tree = inputTree.tree()
            for attr in list(tree.__dict__.keys()):
                if attr.startswith("b_out_"):
                    delattr(tree, attr)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        if not self._is_jetvar_mode():
            for objName in ["Lepton1", "Lepton2", "Lepton3"]:
                for varName in [
                    "pt", "eta", "phi", "mass", "PFIso", "jetPtRatio",
                    "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "jetBTagDeepFlavor"
                ]:
                    self.out.branch("%s_%s" % (objName, varName), "F")
            for objName in ["ptsortLepton1", "ptsortLepton2", "ptsortLepton3"]:
                for varName in ["pt", "eta", "phi"]:
                    self.out.branch("%s_%s" % (objName, varName), "F")
            for objName in ["Z", "LeadingElectron", "LeadingMuon", "LeadingTau"]:
                for varName in ["pt", "eta", "phi", "mass"]:
                    self.out.branch("%s_%s" % (objName, varName), "F")
            self.out.branch("MET_pt", "F")
            self.out.branch("MET_phi", "F")
            self.out.branch("Lepton1_pdgId", "I")
            self.out.branch("Lepton2_pdgId", "I")
            self.out.branch("Lepton3_pdgId", "I")
            self.out.branch("Lepton1_tauGenPartFlav", "I")
            self.out.branch("Lepton2_tauGenPartFlav", "I")
            self.out.branch("Lepton3_tauGenPartFlav", "I")
            self.out.branch("Lepton1_tauDecayMode", "I")
            self.out.branch("Lepton2_tauDecayMode", "I")
            self.out.branch("Lepton3_tauDecayMode", "I")
            self.out.branch("Lepton1_tauProng", "I")
            self.out.branch("Lepton2_tauProng", "I")
            self.out.branch("Lepton3_tauProng", "I")
            self.out.branch("Lepton1_isTight", "O")
            self.out.branch("Lepton2_isTight", "O")
            self.out.branch("Lepton3_isTight", "O")
            self.out.branch("Lepton1_isFake", "O")
            self.out.branch("Lepton2_isFake", "O")
            self.out.branch("Lepton3_isFake", "O")
            self.out.branch("Lepton1_mvaTOP", "F")
            self.out.branch("Lepton2_mvaTOP", "F")
            self.out.branch("Lepton3_mvaTOP", "F")
            self.out.branch("LeadingLepton_pt", "F")
            self.out.branch("TriLepton_mass", "F")
            self.out.branch("TriLepton_pt", "F")
            self.out.branch("TriLepton_WleptonZdPhi", "F")
            self.out.branch("TriLepton_WleptonZdR", "F")
            self.out.branch("Z_charge", "I")
            self.out.branch("nVetoLepton", "i")
            self.out.branch("nGoodLepton", "i")
            self.out.branch("nNPLepton", "i")
            self.out.branch("GoodLeptonCode", "I")
            self.out.branch("nGoodElectron", "i")
            self.out.branch("nGoodMuon", "i")
            self.out.branch("nGoodTau", "i")
            self.out.branch("nNPTau", "i")
            self.out.branch("GoodJet_index", "i", lenVar="nGoodJet")
            for varName in ["pt", "eta", "phi", "mass", "DeepFlavB"]:
                self.out.branch("GoodJet_%s" % varName, "F", lenVar="nGoodJet")
            if self.use_hadronFlavour:
                self.out.branch("GoodJet_hadronFlavour", "I", lenVar="nGoodJet")
            for objName in ["ptsortGoodJet1", "ptsortGoodJet2", "ptsortGoodJet3"]:
                for varName in ["pt", "eta", "phi"]:
                    self.out.branch("%s_%s" % (objName, varName), "F")
            self.out.branch("nBjet", "i")
            self.out.branch("W_MT", "F")
            self.out.branch("nCjet", "i")
            self.out.branch("nUDSjet", "i")
            self.out.branch("GoodJet_DeepFlavCvL", "F", lenVar="nGoodJet")
            self.out.branch("GoodJet_DeepFlavCvB", "F", lenVar="nGoodJet")
            self.out.branch("GoodJet_DeepFlavQG", "F", lenVar="nGoodJet")
        else:
            suffix = self._jet_suffix()
            lenVar = "nGoodJet_%s" % suffix
            self.out.branch("MET_pt_%s" % suffix, "F")
            self.out.branch("MET_phi_%s" % suffix, "F")
            self.out.branch("GoodJet_index_%s" % suffix, "i", lenVar=lenVar)
            for varName in ["pt", "eta", "phi", "mass", "DeepFlavB"]:
                self.out.branch("GoodJet_%s_%s" % (varName, suffix), "F", lenVar=lenVar)
            self.out.branch("GoodJet_hadronFlavour_%s" % suffix, "I", lenVar=lenVar)
            for objName in ["ptsortGoodJet1", "ptsortGoodJet2", "ptsortGoodJet3"]:
                for varName in ["pt", "eta", "phi"]:
                    self.out.branch("%s_%s_%s" % (objName, varName, suffix), "F")
            self.out.branch("nBjet_%s" % suffix, "i")
            self.out.branch("W_MT_%s" % suffix, "F")
            self.out.branch("nCjet_%s" % suffix, "i")
            self.out.branch("nUDSjet_%s" % suffix, "i")
            self.out.branch("GoodJet_DeepFlavCvL_%s" % suffix, "F", lenVar=lenVar)
            self.out.branch("GoodJet_DeepFlavCvB_%s" % suffix, "F", lenVar=lenVar)
            self.out.branch("GoodJet_DeepFlavQG_%s" % suffix, "F", lenVar=lenVar)

        self.initReaders(inputTree)
    def initReaders(self, tree):
        # NanoAODv15 stores the Run-3 prompt lepton discriminator as
        # Electron_promptMVA / Muon_promptMVA.  Older private NanoAODs used
        # Electron_mvaTOP / Muon_mvaTOP.  Keep the internal/output name
        # "mvaTOP" for backward compatibility, but read promptMVA first.
        self.b_MET_pt = self._value_reader_first(tree, ["PuppiMET_pt", "MET_T1_pt", "MET_pt", "PFMET_pt"])
        self.b_MET_phi = self._value_reader_first(tree, ["PuppiMET_phi", "MET_T1_phi", "MET_phi", "PFMET_phi"])

        objName = "Electron"
        setattr(self, "b_n%s" % objName, tree.valueReader("n%s" % objName))
        for varName in [
            "pt", "eta", "isPFcand", "lostHits", "dxy", "dz", "sip3d",
            "miniPFRelIso_all", "sieie", "hoe", "eInvMinusPInv", "convVeto",
            "phi", "mass", "charge", "jetRelIso", "pfRelIso03_all",
            "superclusterEta", "deltaEtaSC"
        ]:
            setattr(self, "b_%s_%s" % (objName, varName), tree.arrayReader("%s_%s" % (objName, varName)))
        self.b_Electron_cutBased = self._array_reader_first(tree, ["Electron_cutBased"], required=False, default=4)
        self.b_Electron_tightCharge = self._array_reader_first(tree, ["Electron_tightCharge"], required=False, default=2)
        self.b_Electron_mvaTOP = self._array_reader_first(tree, ["Electron_promptMVA", "Electron_mvaTOP"])
        self.b_Electron_jetIdx = self._array_reader_first(tree, ["Electron_jetIdx"], required=False, default=-1)
        self.b_Electron_jetBTagDeepFlavor = self._array_reader_first(
            tree, ["Electron_jetDF", "Electron_jetBTagDeepFlavor"], required=False, default=0.0
        )
        self.b_Electron_mvaFall17V2noIso_WPL = self._array_reader_first(
            tree, ["Electron_mvaFall17V2noIso_WPL", "Electron_mvaNoIso_WP90"], required=False, default=True
        )

        objName = "Muon"
        setattr(self, "b_n%s" % objName, tree.valueReader("n%s" % objName))
        for varName in [
            "pt", "eta", "isPFcand", "isTracker", "isGlobal", "dxy", "dz", "sip3d",
            "miniPFRelIso_all", "mediumId", "phi", "mass", "charge", "pfRelIso04_all"
        ]:
            setattr(self, "b_%s_%s" % (objName, varName), tree.arrayReader("%s_%s" % (objName, varName)))
        self.b_Muon_looseId = self._array_reader_first(tree, ["Muon_looseId"], required=False, default=True)
        self.b_Muon_mvaTOP = self._array_reader_first(tree, ["Muon_promptMVA", "Muon_mvaTOP"])
        self.b_Muon_jetIdx = self._array_reader_first(tree, ["Muon_jetIdx"], required=False, default=-1)
        self.b_Muon_jetRelIso = self._array_reader_first(tree, ["Muon_jetRelIso"], required=False, default=-1.0)
        self.b_Muon_jetBTagDeepFlavor = self._array_reader_first(
            tree, ["Muon_jetDF", "Muon_jetBTagDeepFlavor"], required=False, default=0.0
        )

        objName = "Tau"
        setattr(self, "b_n%s" % objName, tree.valueReader("n%s" % objName))
        for varName in ["pt", "eta", "phi", "mass", "charge", "dz", "decayMode", "rawIso", "jetIdx"]:
            # Tau_decayMode is UChar_t in NanoAODv15; convert it with
            # _uchar_to_int() whenever it is used numerically.
            setattr(self, "b_%s_%s" % (objName, varName), tree.arrayReader("%s_%s" % (objName, varName)))
        # NanoAODv15 2024 stores DeepTau2018v2p5 branches.
        self.b_Tau_idDeepTau2018v2p5VSjet = tree.arrayReader("Tau_idDeepTau2018v2p5VSjet")
        self.b_Tau_idDeepTau2018v2p5VSe = tree.arrayReader("Tau_idDeepTau2018v2p5VSe")
        self.b_Tau_idDeepTau2018v2p5VSmu = tree.arrayReader("Tau_idDeepTau2018v2p5VSmu")
        raw_tree = self._raw_tree(tree)
        self.has_Tau_genPartFlav = bool(hasattr(raw_tree, "GetBranch") and raw_tree.GetBranch("Tau_genPartFlav"))
        if self.has_Tau_genPartFlav:
            self.b_Tau_genPartFlav = tree.arrayReader("Tau_genPartFlav")

        objName = "Jet"
        setattr(self, "b_n%s" % objName, tree.valueReader("n%s" % objName))

        # Central NanoAODv15 2024 AK4PUPPI jet four-vector branches.
        for varName in ["pt", "eta", "phi", "mass"]:
            setattr(self, "b_%s_%s" % (objName, varName), tree.arrayReader("%s_%s" % (objName, varName)))

        # Variables needed to emulate tight jet ID for NanoAODv15 PUPPI jets.
        # They are present in both the 2024 data and MC documentation.
        for varName in [
            "chHEF", "neHEF", "chEmEF", "neEmEF", "muEF",
            "chMultiplicity", "neMultiplicity", "nConstituents",
        ]:
            setattr(
                self,
                "b_%s_%s" % (objName, varName),
                self._array_reader_first(tree, ["%s_%s" % (objName, varName)], required=False, default=0.0),
            )

        # 2024 NanoAODv15 stores the official AK4 BTV discriminator as UParTAK4.
        # Keep the output branch names as GoodJet_DeepFlav* for downstream
        # compatibility, but read the NanoAODv15 UParTAK4 branch names here.
        self.b_Jet_btagDeepFlavB = tree.arrayReader("Jet_btagUParTAK4B")
        self.b_Jet_btagDeepFlavCvL = self._array_reader_first(
            tree, ["Jet_btagUParTAK4CvL"], required=False, default=-999.0
        )
        self.b_Jet_btagDeepFlavCvB = self._array_reader_first(
            tree, ["Jet_btagUParTAK4CvB"], required=False, default=-999.0
        )
        self.b_Jet_btagDeepFlavQG = self._array_reader_first(
            tree, ["Jet_btagUParTAK4QvG"], required=False, default=-999.0
        )

        if self.use_hadronFlavour or self._is_jetvar_mode():
            self.b_Jet_hadronFlavour = self._array_reader_first(
                tree, ["Jet_hadronFlavour"], required=False, default=0
            )

        if self._is_jetvar_mode():
            self.b_Jet_Uncert = self._array_reader_first(
                tree, ["Jet_jecUncert%s" % self.jetVar], required=False, default=0.0
            )

        # Lepton-associated jet four-vectors.  NanoAODv15 no longer provides
        # Electron_jetPt/Muon_jetPt style branches; it provides jetIdx instead.
        if self._has_branch(tree, "Electron_jetPt"):
            self.b_Electron_jetPt = tree.arrayReader("Electron_jetPt")
            self.b_Electron_jetEta = tree.arrayReader("Electron_jetEta")
            self.b_Electron_jetPhi = tree.arrayReader("Electron_jetPhi")
            self.b_Electron_jetMass = tree.arrayReader("Electron_jetMass")
        else:
            self.b_Electron_jetPt = _JetComponentFromIndexReader(self.b_Electron_jetIdx, self.b_Jet_pt, 0.0)
            self.b_Electron_jetEta = _JetComponentFromIndexReader(self.b_Electron_jetIdx, self.b_Jet_eta, 0.0)
            self.b_Electron_jetPhi = _JetComponentFromIndexReader(self.b_Electron_jetIdx, self.b_Jet_phi, 0.0)
            self.b_Electron_jetMass = _JetComponentFromIndexReader(self.b_Electron_jetIdx, self.b_Jet_mass, 0.0)

        if self._has_branch(tree, "Muon_jetPt"):
            self.b_Muon_jetPt = tree.arrayReader("Muon_jetPt")
            self.b_Muon_jetEta = tree.arrayReader("Muon_jetEta")
            self.b_Muon_jetPhi = tree.arrayReader("Muon_jetPhi")
            self.b_Muon_jetMass = tree.arrayReader("Muon_jetMass")
        else:
            self.b_Muon_jetPt = _JetComponentFromIndexReader(self.b_Muon_jetIdx, self.b_Jet_pt, 0.0)
            self.b_Muon_jetEta = _JetComponentFromIndexReader(self.b_Muon_jetIdx, self.b_Jet_eta, 0.0)
            self.b_Muon_jetPhi = _JetComponentFromIndexReader(self.b_Muon_jetIdx, self.b_Jet_phi, 0.0)
            self.b_Muon_jetMass = _JetComponentFromIndexReader(self.b_Muon_jetIdx, self.b_Jet_mass, 0.0)

        if self._has_branch(tree, "Electron_jetPtRatio"):
            self.b_Electron_jetPtRatio = tree.arrayReader("Electron_jetPtRatio")
        else:
            self.b_Electron_jetPtRatio = _JetPtRatioFromRelIsoReader(self.b_Electron_jetRelIso)

        if self._has_branch(tree, "Muon_jetPtRatio"):
            self.b_Muon_jetPtRatio = tree.arrayReader("Muon_jetPtRatio")
        else:
            self.b_Muon_jetPtRatio = _JetPtRatioFromRelIsoReader(self.b_Muon_jetRelIso)

        self._ttreereaderversion = tree._ttreereaderversion
    def setCjetCuts(self, cvl, cvb):
        self.minCjetCvLDiscr_ = float(cvl)
        self.minCjetCvBDiscr_ = float(cvb)

    def setQjetCut(self, qg):
        self.minQjetQGDiscr_ = float(qg)

    def _is_taufr_mode(self):
        return bool(self.mode) and self.mode.startswith("FR")

    def _candidate_pt(self, flavor, idx):
        if flavor == "e":
            return float(self.b_Electron_pt.At(idx))
        if flavor == "mu":
            return float(self.b_Muon_pt.At(idx))
        if flavor == "tau":
            return float(self.b_Tau_pt.At(idx))
        return -999.0

    def _candidate_is_fake(self, flavor, idx, electronIdxs, muonIdxs):
        if flavor == "e":
            return bool(self.isOnlyFakeElectron(idx))
        if flavor == "mu":
            return bool(self.isOnlyFakeMuon(idx))
        if flavor == "tau":
            return bool(self.isOnlyFakeTau(idx, electronIdxs, muonIdxs))
        return False

    def _iter_combinations(self, idxs, n):
        if n <= 0:
            yield ()
            return
        for combo in itertools.combinations(idxs, n):
            yield combo

    def _select_npl_indices(self, mode, npElectronIdxs, npMuonIdxs, npTauIdxs):
        req_e, req_m, req_t = self._required_lepton_multiplicity(mode)
        if len(npElectronIdxs) < req_e or len(npMuonIdxs) < req_m or len(npTauIdxs) < req_t:
            return None

        best = None
        best_score = None

        for ecombo in self._iter_combinations(npElectronIdxs, req_e):
            ecombo = list(ecombo)
            for mucombo in self._iter_combinations(npMuonIdxs, req_m):
                mucombo = list(mucombo)
                for taucombo in self._iter_combinations(npTauIdxs, req_t):
                    taucombo = list(taucombo)

                    fake_count = 0
                    fake_count += sum(1 for i in ecombo if self._candidate_is_fake("e", i, ecombo, mucombo))
                    fake_count += sum(1 for i in mucombo if self._candidate_is_fake("mu", i, ecombo, mucombo))
                    fake_count += sum(1 for i in taucombo if self._candidate_is_fake("tau", i, ecombo, mucombo))
                    if fake_count < 1:
                        continue

                    pts = sorted(
                        [self._candidate_pt("e", i) for i in ecombo] +
                        [self._candidate_pt("mu", i) for i in mucombo] +
                        [self._candidate_pt("tau", i) for i in taucombo],
                        reverse=True
                    )
                    # Prefer the highest-pT valid combination; use fewer replacements (larger pt sum)
                    score = tuple(pts + [fake_count])
                    if best_score is None or score > best_score:
                        best_score = score
                        best = (ecombo, mucombo, taucombo)

        return best

    def _is_jetvar_mode(self):
        return self.jetVar is not None and self.jetSyst in ("Up", "Down")

    def _jet_shift_sign(self):
        if self.jetSyst == "Up":
            return 1.0
        if self.jetSyst == "Down":
            return -1.0
        return 0.0

    def _jet_suffix(self):
        if not self._is_jetvar_mode():
            return ""
        return "%s%s" % (self.jetVar, self.jetSyst)

    def _is_btaggable_jet_for_tag_count(self, i):
        # UParTAK4 b-tagging WPs/SFs are provided for AK4 jets in the tracker
        # b-tagging acceptance.  AN2025/031 uses the loose WP for jets with
        # pT > 20 GeV and |eta| < 2.5.  This threshold is lower than the
        # analysis GoodJet pT threshold (30 GeV), so nBjet is counted from the
        # raw selected-jet collection, not only from GoodJets.
        pt = float(self.b_Jet_pt.At(i))
        if self._is_jetvar_mode():
            pt = pt * (1.0 + self._jet_shift_sign() * float(self.b_Jet_Uncert.At(i)))
        return (pt > 20.0) and (abs(float(self.b_Jet_eta.At(i))) < self.maxBTagJetAbsEta_)

    @staticmethod
    def _p4_from_arrays(pt_arr, eta_arr, phi_arr, mass_arr, idx):
        p4 = ROOT.TLorentzVector()
        p4.SetPtEtaPhiM(pt_arr.At(idx), eta_arr.At(idx), phi_arr.At(idx), mass_arr.At(idx))
        return p4

    @staticmethod
    def _uchar_to_int(x):
        """Convert ROOT UChar_t/Char_t array values to Python int.

        In NanoAODTools with Python 3, a UChar_t branch can be returned as a
        one-character string, e.g. Tau_decayMode == "\x0b".  Calling
        int("\x0b") raises ValueError, so all UChar_t-like branches must
        pass through this converter before numerical comparisons or branch
        filling.
        """
        if isinstance(x, str):
            if len(x) == 0:
                return 0
            return ord(x[0])
        if isinstance(x, bytes):
            if len(x) == 0:
                return 0
            return x[0]
        return int(x)

    @staticmethod
    def _copy_p4_from_lv(lv):
        return [float(lv.Pt()), float(lv.Eta()), float(lv.Phi()), float(lv.M())]

    @staticmethod
    def _copy_p4_from_arrays(pt_arr, eta_arr, phi_arr, mass_arr, idx):
        return [
            float(pt_arr.At(idx)),
            float(eta_arr.At(idx)),
            float(phi_arr.At(idx)),
            float(mass_arr.At(idx)),
        ]

    @staticmethod
    def _compute_mt(lepP4, met_pt, met_phi):
        met_px = met_pt * math.cos(met_phi)
        met_py = met_pt * math.sin(met_phi)
        pt = lepP4.Pt() + met_pt
        px = lepP4.Px() + met_px
        py = lepP4.Py() + met_py
        return math.sqrt(max(0.0, pt * pt - px * px - py * py))

    @staticmethod
    def _set_tree_attr(tree, name, value):
        setattr(tree, "b_out_%s" % name, value)

    @staticmethod
    def _get_tree_attr(tree, name, default=None):
        attr_name = "b_out_%s" % name
        if hasattr(tree, attr_name):
            return getattr(tree, attr_name)
        return default

    def _fill_branch_and_cache(self, tree, name, value):
        self._set_tree_attr(tree, name, value)
        self.out.fillBranch(name, value)

    def _get_central_lepton_p4(self, tree, idx):
        p4 = ROOT.TLorentzVector()
        p4.SetPtEtaPhiM(
            float(self._get_tree_attr(tree, f"Lepton{idx}_pt", 0.0)),
            float(self._get_tree_attr(tree, f"Lepton{idx}_eta", 0.0)),
            float(self._get_tree_attr(tree, f"Lepton{idx}_phi", 0.0)),
            float(self._get_tree_attr(tree, f"Lepton{idx}_mass", 0.0)),
        )
        return p4

    @staticmethod
    def _get_ptsort_goodjet(goodjet_p4, nGoodJet, which, component):
        if nGoodJet >= which:
            return goodjet_p4[component][which - 1]
        return -999.0

    @staticmethod
    def _tau_prong_group_from_decaymode(decay_mode):
        dm = FCNCTriLepton._uchar_to_int(decay_mode)
        if dm in (0, 1, 2):
            return 1
        if dm in (10, 11):
            return 23
        if dm >= 10:
            return 23
        return 0

    @staticmethod
    def _required_lepton_multiplicity(mode):
        """
        mapping = {
            "MuElEl": (1, 2, 0),
            "ElMuMu": (2, 1, 0),
            "ElElEl": (3, 0, 0),
            "MuMuMu": (0, 3, 0),
            "TaTaTa": (0, 0, 3),
            "ElTaTa": (1, 0, 2),
            "MuTaTa": (0, 1, 2),
            "TaElEl": (2, 0, 1),
            "TaMuMu": (0, 2, 1),
            "NPLMuElEl": (1, 2, 0),
            "NPLElMuMu": (2, 1, 0),
            "NPLElElEl": (3, 0, 0),
            "NPLMuMuMu": (0, 3, 0),
            "NPLTaTaTa": (0, 0, 3),
            "NPLElTaTa": (1, 0, 2),
            "NPLMuTaTa": (0, 1, 2),
            "NPLTaElEl": (2, 0, 1),
            "NPLTaMuMu": (0, 2, 1),
        }
        """
        mapping = {
            "MuElEl": (2, 1, 0),
            "ElMuMu": (1, 2, 0),
            "ElElEl": (3, 0, 0),
            "MuMuMu": (0, 3, 0),
            "TaTaTa": (0, 0, 3),
            "ElTaTa": (1, 0, 2),
            "MuTaTa": (0, 1, 2),
            "TaElEl": (2, 0, 1),
            "TaMuMu": (0, 2, 1),
            "NPLMuElEl": (2, 1, 0),
            "NPLElMuMu": (1, 2, 0),
            "NPLElElEl": (3, 0, 0),
            "NPLMuMuMu": (0, 3, 0),
            "NPLTaTaTa": (0, 0, 3),
            "NPLElTaTa": (1, 0, 2),
            "NPLMuTaTa": (0, 1, 2),
            "NPLTaElEl": (2, 0, 1),
            "NPLTaMuMu": (0, 2, 1),
        }
        return mapping.get(mode, (0, 0, 0))

    def _get_lepton_aux_state(self, idx):
        return {
            "pdgId": getattr(self, f"out_Lepton{idx}_pdgId"),
            "tauDecayMode": getattr(self, f"out_Lepton{idx}_tauDecayMode"),
            "tauProng": getattr(self, f"out_Lepton{idx}_tauProng"),
            "tauGenPartFlav": getattr(self, f"out_Lepton{idx}_tauGenPartFlav"),
            "PFIso": getattr(self, f"out_Lepton{idx}_PFIso"),
            "jetPtRatio": getattr(self, f"out_Lepton{idx}_jetPtRatio"),
            "isTight": getattr(self, f"out_Lepton{idx}_isTight"),
            "isFake": getattr(self, f"out_Lepton{idx}_isFake"),
            "jetBTagDeepFlavor": getattr(self, f"out_Lepton{idx}_jetBTagDeepFlavor"),
            "mvaTOP": getattr(self, f"out_Lepton{idx}_mvaTOP"),
            "Jet_p4": list(getattr(self, f"out_Lepton{idx}_Jet_p4")),
        }

    def _set_lepton_aux_state(self, idx, state):
        setattr(self, f"out_Lepton{idx}_pdgId", int(state["pdgId"]))
        setattr(self, f"out_Lepton{idx}_tauDecayMode", int(state["tauDecayMode"]))
        setattr(self, f"out_Lepton{idx}_tauProng", int(state["tauProng"]))
        setattr(self, f"out_Lepton{idx}_tauGenPartFlav", int(state["tauGenPartFlav"]))
        setattr(self, f"out_Lepton{idx}_PFIso", float(state["PFIso"]))
        setattr(self, f"out_Lepton{idx}_jetPtRatio", float(state["jetPtRatio"]))
        setattr(self, f"out_Lepton{idx}_isTight", bool(state["isTight"]))
        setattr(self, f"out_Lepton{idx}_isFake", bool(state["isFake"]))
        setattr(self, f"out_Lepton{idx}_jetBTagDeepFlavor", float(state["jetBTagDeepFlavor"]))
        setattr(self, f"out_Lepton{idx}_mvaTOP", float(state["mvaTOP"]))
        setattr(self, f"out_Lepton{idx}_Jet_p4", list(state["Jet_p4"]))

    def _permute_lepton_aux_states(self, order):
        old_states = [self._get_lepton_aux_state(i) for i in (1, 2, 3)]
        for new_idx, old_idx in enumerate(order, start=1):
            self._set_lepton_aux_state(new_idx, old_states[old_idx])

    @staticmethod
    def _sort_desc_by_pt(idxs, pt_reader):
        idxs.sort(key=lambda i: pt_reader.At(i), reverse=True)

    def _get_corrected_met(self):
        met_px = self.out_MET_pt * math.cos(self.out_MET_phi)
        met_py = self.out_MET_pt * math.sin(self.out_MET_phi)
        shift_sign = self._jet_shift_sign()
        for pt, phi, uncert in zip(self.out_GoodJet_p4[0], self.out_GoodJet_p4[2], self.out_GoodJet_Uncert):
            met_px += (-shift_sign) * uncert * pt * math.cos(phi)
            met_py += (-shift_sign) * uncert * pt * math.sin(phi)
        return float(math.hypot(met_px, met_py)), float(math.atan2(met_py, met_px))

    def _get_shifted_goodjet_pt(self):
        shift_sign = self._jet_shift_sign()
        return [float(pt * (1.0 + shift_sign * uncert)) for pt, uncert in zip(self.out_GoodJet_p4[0], self.out_GoodJet_Uncert)]

    def _get_shifted_goodjet_mass(self):
        shift_sign = self._jet_shift_sign()
        return [float(mass * (1.0 + shift_sign * uncert)) for mass, uncert in zip(self.out_GoodJet_p4[3], self.out_GoodJet_Uncert)]

    def _analyze_jetvar_only(self, event):
        tree = event._tree

        if self._get_tree_attr(tree, "GoodLeptonCode", None) is None:
            raise RuntimeError(
                'jetVar mode requires central fcncTriLepton_allmodes_unified to run first in the same module chain.'
            )

        self.out_GoodLeptonCode = int(self._get_tree_attr(tree, "GoodLeptonCode", 0))
        if abs(self.out_GoodLeptonCode) < 100:
            return False

        lepton1P4 = self._get_central_lepton_p4(tree, 1)
        lepton2P4 = self._get_central_lepton_p4(tree, 2)
        lepton3P4 = self._get_central_lepton_p4(tree, 3)

        self.out_MET_pt = float(self._get_tree_attr(tree, "MET_pt", self._event_value_first(event, ["PuppiMET_pt", "MET_T1_pt", "MET_pt", "PFMET_pt"])))
        self.out_MET_phi = float(self._get_tree_attr(tree, "MET_phi", self._event_value_first(event, ["PuppiMET_phi", "MET_T1_phi", "MET_phi", "PFMET_phi"])))

        jetIdxs = []
        for i in range(self.b_Jet_pt.GetSize()):
            if not self.isGoodJet(i):
                continue
            jetP4 = self._p4_from_arrays(self.b_Jet_pt, self.b_Jet_eta, self.b_Jet_phi, self.b_Jet_mass, i)
            if lepton1P4.Pt() > 0 and lepton1P4.DeltaR(jetP4) < 0.4:
                continue
            if lepton2P4.Pt() > 0 and lepton2P4.DeltaR(jetP4) < 0.4:
                continue
            if lepton3P4.Pt() > 0 and lepton3P4.DeltaR(jetP4) < 0.4:
                continue
            jetIdxs.append(i)

            isTaggableJet = self._is_btaggable_jet_for_tag_count(i)
            isBTagged = isTaggableJet and self.b_Jet_btagDeepFlavB.At(i) > self.minBjetBDiscr_
            isCTagged = isTaggableJet and (
                self.b_Jet_btagDeepFlavCvL.At(i) > self.minCjetCvLDiscr_ and
                self.b_Jet_btagDeepFlavCvB.At(i) > self.minCjetCvBDiscr_
            )
            isQTagged = isTaggableJet and self.b_Jet_btagDeepFlavQG.At(i) > self.minQjetQGDiscr_

            if isBTagged:
                self.out_nBjet += 1
            if isCTagged and not isBTagged:
                self.out_nCjet += 1
            if isQTagged and not isBTagged and not isCTagged:
                self.out_nUDSjet += 1

        self.out_nGoodJet = len(jetIdxs)
        jetIdxs.sort(key=lambda i: self.b_Jet_pt.At(i), reverse=True)
        for kk in jetIdxs:
            self.out_GoodJet_p4[0].append(float(self.b_Jet_pt.At(kk)))
            self.out_GoodJet_p4[1].append(float(self.b_Jet_eta.At(kk)))
            self.out_GoodJet_p4[2].append(float(self.b_Jet_phi.At(kk)))
            self.out_GoodJet_p4[3].append(float(self.b_Jet_mass.At(kk)))
            self.out_GoodJet_DeepFlavB.append(float(self.b_Jet_btagDeepFlavB.At(kk)))
            self.out_GoodJet_index.append(int(kk))
            self.out_GoodJet_DeepFlavCvL.append(float(self.b_Jet_btagDeepFlavCvL.At(kk)))
            self.out_GoodJet_DeepFlavCvB.append(float(self.b_Jet_btagDeepFlavCvB.At(kk)))
            self.out_GoodJet_DeepFlavQG.append(float(self.b_Jet_btagDeepFlavQG.At(kk)))
            self.out_GoodJet_hadronFlavour.append(self._uchar_to_int(self.b_Jet_hadronFlavour.At(kk)))
            self.out_GoodJet_Uncert.append(float(self.b_Jet_Uncert.At(kk)))

        suffix = self._jet_suffix()
        corrected_met_pt, corrected_met_phi = self._get_corrected_met()
        shifted_w_mt = self._compute_mt(lepton1P4, corrected_met_pt, corrected_met_phi)

        ptsort_goodjet_map = {
            "ptsortGoodJet1_pt_%s" % suffix: self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 1, 0),
            "ptsortGoodJet1_eta_%s" % suffix: self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 1, 1),
            "ptsortGoodJet1_phi_%s" % suffix: self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 1, 2),
            "ptsortGoodJet2_pt_%s" % suffix: self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 2, 0),
            "ptsortGoodJet2_eta_%s" % suffix: self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 2, 1),
            "ptsortGoodJet2_phi_%s" % suffix: self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 2, 2),
            "ptsortGoodJet3_pt_%s" % suffix: self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 3, 0),
            "ptsortGoodJet3_eta_%s" % suffix: self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 3, 1),
            "ptsortGoodJet3_phi_%s" % suffix: self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 3, 2),
        }
        for name, value in ptsort_goodjet_map.items():
            self._fill_branch_and_cache(tree, name, value)

        other_map = {
            "GoodJet_pt_%s" % suffix: self._get_shifted_goodjet_pt(),
            "GoodJet_eta_%s" % suffix: self.out_GoodJet_p4[1],
            "GoodJet_phi_%s" % suffix: self.out_GoodJet_p4[2],
            "GoodJet_mass_%s" % suffix: self._get_shifted_goodjet_mass(),
            "W_MT_%s" % suffix: shifted_w_mt,
            "GoodJet_index_%s" % suffix: self.out_GoodJet_index,
            "GoodJet_DeepFlavB_%s" % suffix: self.out_GoodJet_DeepFlavB,
            "GoodJet_hadronFlavour_%s" % suffix: self.out_GoodJet_hadronFlavour,
            "nBjet_%s" % suffix: self.out_nBjet,
            "nCjet_%s" % suffix: self.out_nCjet,
            "nUDSjet_%s" % suffix: self.out_nUDSjet,
            "GoodJet_DeepFlavCvL_%s" % suffix: self.out_GoodJet_DeepFlavCvL,
            "GoodJet_DeepFlavCvB_%s" % suffix: self.out_GoodJet_DeepFlavCvB,
            "GoodJet_DeepFlavQG_%s" % suffix: self.out_GoodJet_DeepFlavQG,
            "MET_pt_%s" % suffix: corrected_met_pt,
            "MET_phi_%s" % suffix: corrected_met_phi,
            "nGoodJet_%s" % suffix: self.out_nGoodJet,
        }
        for name, value in other_map.items():
            self._fill_branch_and_cache(tree, name, value)

        return True

    def _reset_outputs(self):
        self.out_LeadingMuon_p4 = [0.0, 0.0, 0.0, 0.0]
        self.out_LeadingElectron_p4 = [0.0, 0.0, 0.0, 0.0]
        self.out_LeadingTau_p4 = [0.0, 0.0, 0.0, 0.0]

        self.out_Lepton1_p4 = [0.0, 0.0, 0.0, 0.0]
        self.out_Lepton2_p4 = [0.0, 0.0, 0.0, 0.0]
        self.out_Lepton3_p4 = [0.0, 0.0, 0.0, 0.0]

        self.out_Lepton1_Jet_p4 = [0.0, 0.0, 0.0, 0.0]
        self.out_Lepton2_Jet_p4 = [0.0, 0.0, 0.0, 0.0]
        self.out_Lepton3_Jet_p4 = [0.0, 0.0, 0.0, 0.0]

        self.out_Lepton1_pdgId = 0
        self.out_Lepton2_pdgId = 0
        self.out_Lepton3_pdgId = 0
        self.out_Lepton1_tauDecayMode = -1
        self.out_Lepton2_tauDecayMode = -1
        self.out_Lepton3_tauDecayMode = -1
        self.out_Lepton1_tauProng = 0
        self.out_Lepton2_tauProng = 0
        self.out_Lepton3_tauProng = 0
        self.out_Lepton1_tauGenPartFlav = -1
        self.out_Lepton2_tauGenPartFlav = -1
        self.out_Lepton3_tauGenPartFlav = -1

        self.out_Lepton1_PFIso = 0.0
        self.out_Lepton2_PFIso = 0.0
        self.out_Lepton3_PFIso = 0.0

        self.out_Lepton1_jetPtRatio = 0.0
        self.out_Lepton2_jetPtRatio = 0.0
        self.out_Lepton3_jetPtRatio = 0.0

        self.out_Lepton1_isTight = False
        self.out_Lepton2_isTight = False
        self.out_Lepton3_isTight = False

        self.out_Lepton1_isFake = False
        self.out_Lepton2_isFake = False
        self.out_Lepton3_isFake = False

        self.out_Lepton1_jetBTagDeepFlavor = 0.0
        self.out_Lepton2_jetBTagDeepFlavor = 0.0
        self.out_Lepton3_jetBTagDeepFlavor = 0.0

        self.out_Lepton1_mvaTOP = 0.0
        self.out_Lepton2_mvaTOP = 0.0
        self.out_Lepton3_mvaTOP = 0.0

        self.out_ptsortLepton1_pt = 0.0
        self.out_ptsortLepton1_eta = 0.0
        self.out_ptsortLepton1_phi = 0.0
        self.out_ptsortLepton2_pt = 0.0
        self.out_ptsortLepton2_eta = 0.0
        self.out_ptsortLepton2_phi = 0.0
        self.out_ptsortLepton3_pt = 0.0
        self.out_ptsortLepton3_eta = 0.0
        self.out_ptsortLepton3_phi = 0.0

        self.out_TriLepton_mass = 0.0
        self.out_TriLepton_pt = 0.0
        self.out_TriLepton_WleptonZdPhi = 0.0
        self.out_TriLepton_WleptonZdR = 0.0

        self.out_Z_p4 = [0.0, 0.0, 0.0, 0.0]
        self.out_Z_charge = 0

        self.out_MET_pt = 0.0
        self.out_MET_phi = 0.0
        self.out_W_MT = 0.0

        self.out_GoodLeptonCode = 0
        self.out_nGoodElectron = 0
        self.out_nGoodMuon = 0
        self.out_nGoodTau = 0
        self.out_nNPElectron = 0
        self.out_nNPMuon = 0
        self.out_nNPTau = 0
        self.out_nVetoElectron = 0
        self.out_nVetoMuon = 0

        self.out_nGoodJet = 0
        self.out_nBjet = 0
        self.out_nCjet = 0
        self.out_nUDSjet = 0

        self.out_GoodJet_p4 = [[], [], [], []]
        self.out_GoodJet_DeepFlavB = []
        self.out_GoodJet_index = []
        self.out_GoodJet_DeepFlavCvL = []
        self.out_GoodJet_DeepFlavCvB = []
        self.out_GoodJet_DeepFlavQG = []
        self.out_GoodJet_hadronFlavour = []
        self.out_GoodJet_Uncert = []

    # ------------------------------------------------------------------
    # Ported object selections
    # ------------------------------------------------------------------
    def _passes_muon_veto(self, i):
        if not (self.b_Muon_pt.At(i) > 10):
            return False
        if not (abs(self.b_Muon_eta.At(i)) < 2.4):
            return False
        if not self.b_Muon_looseId.At(i):
            return False
        return True

    def _passes_muon_relax(self, i):
        # AN2025/031 Table 8: relax muon = medium ID + impact parameter
        # cuts + isolation/pT < 0.4 + matched-jet relative isolation < 0.5.
        if not self._passes_muon_veto(i):
            return False
        if not self.b_Muon_mediumId.At(i):
            return False
        if not (abs(self.b_Muon_dxy.At(i)) < 0.05):
            return False
        if not (abs(self.b_Muon_dz.At(i)) < 0.10):
            return False
        if not (self.b_Muon_pfRelIso04_all.At(i) < 0.4):
            return False
        jet_rel_iso = float(self.b_Muon_jetRelIso.At(i))
        if jet_rel_iso >= 0.0 and not (jet_rel_iso < 0.5):
            return False
        return True

    def isGoodMuon(self, i):
        # AN2025/031 Table 8 strict muon: relax muon + promptMVA > 0.7.
        if not self._passes_muon_relax(i):
            return False
        if not (self.b_Muon_mvaTOP.At(i) > self.minMuonPromptMVA_):
            return False
        return True

    def isNPMuon(self, i):
        # Nonprompt application denominator = relax muon.  Tight muons are
        # intentionally included; isOnlyFakeMuon marks relax-not-strict legs.
        return self._passes_muon_relax(i)

    def isOnlyFakeMuon(self, i):
        return self.isNPMuon(i) and (not self.isGoodMuon(i))

    def _passes_electron_veto(self, i):
        if not (self.b_Electron_pt.At(i) > 10):
            return False
        if not (abs(self.b_Electron_eta.At(i)) < 2.5):
            return False
        if not (self._uchar_to_int(self.b_Electron_cutBased.At(i)) >= 1):
            return False
        return True

    def _passes_electron_relax(self, i):
        # AN2025/031 Table 9: relax electron = medium ID+isolation + impact
        # parameter cuts + matched-jet relative isolation < 0.5 + conversion veto.
        if not self._passes_electron_veto(i):
            return False
        if not (self._uchar_to_int(self.b_Electron_cutBased.At(i)) >= 3):
            return False
        if not (abs(self.b_Electron_dxy.At(i)) < 0.05):
            return False
        if not (abs(self.b_Electron_dz.At(i)) < 0.10):
            return False
        jet_rel_iso = float(self.b_Electron_jetRelIso.At(i))
        if jet_rel_iso >= 0.0 and not (jet_rel_iso < 0.5):
            return False
        if not self.b_Electron_convVeto.At(i):
            return False
        return True

    def isGoodElectron(self, i):
        # AN2025/031 Table 9 strict electron: relax electron + triple-charge
        # consistency + promptMVA threshold.  For 2024 NanoAODv15 the threshold
        # is 0.9; for 2022-2023 pass electronPromptMVACut=0.3 to the factory.
        if not self._passes_electron_relax(i):
            return False
        if not (self._uchar_to_int(self.b_Electron_tightCharge.At(i)) >= 2):
            return False
        if not (self.b_Electron_mvaTOP.At(i) > self.minElectronPromptMVA_):
            return False
        return True

    def isNPElectron(self, i):
        # Nonprompt application denominator = relax electron.  Tight electrons
        # are intentionally included; isOnlyFakeElectron marks relax-not-strict legs.
        return self._passes_electron_relax(i)

    def isOnlyFakeElectron(self, i):
        return self.isNPElectron(i) and (not self.isGoodElectron(i))

    @staticmethod
    def _deep_tau_pass(score, bit):
        value = FCNCTriLepton._uchar_to_int(score)
        bit = int(bit)
        # Run-2 NanoAOD DeepTau WPs were commonly treated as bit masks.
        # NanoAODv15 DeepTau2018v2p5 branches are documented as ordered WP
        # integers: 1=VVVLoose, ..., 5=Medium, 6=Tight, ... .  Accept both
        # encodings so tau modes do not break when the same module is run on
        # Run-3 NanoAODv15.
        if 0 <= value <= 8 and bit > 0 and (bit & (bit - 1)) == 0:
            wp_index = int(math.log(bit, 2)) + 1
            return value >= wp_index
        return (value & bit) != 0

    def _tau_jet_idx(self, i):
        idx = int(self.b_Tau_jetIdx.At(i))
        if idx < 0 or idx >= self.b_Jet_pt.GetSize():
            return -1
        return idx

    def _tau_matched_jet_p4(self, i):
        jet_idx = self._tau_jet_idx(i)
        if jet_idx < 0:
            return ROOT.TLorentzVector()
        return self._p4_from_arrays(self.b_Jet_pt, self.b_Jet_eta, self.b_Jet_phi, self.b_Jet_mass, jet_idx)

    def _tau_jet_pt_ratio(self, i):
        jet_idx = self._tau_jet_idx(i)
        if jet_idx < 0:
            return 1.0
        jet_pt = float(self.b_Jet_pt.At(jet_idx))
        if jet_pt <= 0.0:
            return 1.0
        return float(self.b_Tau_pt.At(i) / jet_pt)

    def _tau_matched_jet_btag(self, i):
        jet_idx = self._tau_jet_idx(i)
        if jet_idx < 0:
            return 0.0
        return float(self.b_Jet_btagDeepFlavB.At(jet_idx))

    def _tau_passes_lepton_cleaning(self, i, electronIdxs=None, muonIdxs=None):
        tauP4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
        for ie in (electronIdxs or []):
            eleP4 = self._p4_from_arrays(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_phi, self.b_Electron_mass, ie)
            if tauP4.DeltaR(eleP4) < 0.4:
                return False
        for im in (muonIdxs or []):
            muP4 = self._p4_from_arrays(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_phi, self.b_Muon_mass, im)
            if tauP4.DeltaR(muP4) < 0.4:
                return False
        return True

    def _tau_decaymode(self, i):
        return self._uchar_to_int(self.b_Tau_decayMode.At(i))

    def _passes_tau_baseline(self, i, electronIdxs=None, muonIdxs=None):
        if not (self.b_Tau_pt.At(i) > 20):
            return False
        if not (abs(self.b_Tau_eta.At(i)) < 2.3):
            return False
        if not (abs(self.b_Tau_dz.At(i)) < 0.2):
            return False
        if self._tau_decaymode(i) in (5, 6):
            return False
        if not self._tau_passes_lepton_cleaning(i, electronIdxs, muonIdxs):
            return False
        return True

    def _tau_genpartflav(self, i):
        if hasattr(self, "b_Tau_genPartFlav") and self.b_Tau_genPartFlav is not None:
            return self._uchar_to_int(self.b_Tau_genPartFlav.At(i))
        return -1

    def isGoodTau(self, i, electronIdxs=None, muonIdxs=None):
        if not self._passes_tau_baseline(i, electronIdxs, muonIdxs):
            return False
        if not self._deep_tau_pass(self.b_Tau_idDeepTau2018v2p5VSjet.At(i), 16):
            return False
        if not self._deep_tau_pass(self.b_Tau_idDeepTau2018v2p5VSe.At(i), 2):
            return False
        if not self._deep_tau_pass(self.b_Tau_idDeepTau2018v2p5VSmu.At(i), 1):
            return False
        return True

    def isNPTau(self, i, electronIdxs=None, muonIdxs=None):
        if not self._passes_tau_baseline(i, electronIdxs, muonIdxs):
            return False
        if not self._deep_tau_pass(self.b_Tau_idDeepTau2018v2p5VSjet.At(i), 2):
            return False
        if not self._deep_tau_pass(self.b_Tau_idDeepTau2018v2p5VSe.At(i), 1):
            return False
        if not self._deep_tau_pass(self.b_Tau_idDeepTau2018v2p5VSmu.At(i), 1):
            return False
        return True

    def isOnlyFakeTau(self, i, electronIdxs=None, muonIdxs=None):
        return self.isNPTau(i, electronIdxs, muonIdxs) and (not self.isGoodTau(i, electronIdxs, muonIdxs))

    def _passes_official_ak4puppi_tight_jet_id(self, i):
        if self.jetIDCorr_ is None:
            return None
        try:
            ch_mult = self._uchar_to_int(self.b_Jet_chMultiplicity.At(i))
            ne_mult = self._uchar_to_int(self.b_Jet_neMultiplicity.At(i))
            multiplicity = int(ch_mult + ne_mult)
            value = self.jetIDCorr_.evaluate(
                float(self.b_Jet_eta.At(i)),
                float(self.b_Jet_chHEF.At(i)),
                float(self.b_Jet_neHEF.At(i)),
                float(self.b_Jet_chEmEF.At(i)),
                float(self.b_Jet_neEmEF.At(i)),
                float(self.b_Jet_muEF.At(i)),
                int(ch_mult),
                int(ne_mult),
                int(multiplicity),
            )
            return float(value) > 0.5
        except Exception:
            return None

    def _passes_tight_jet_id_from_composition(self, i):
        # Fallback for NanoAODv15 PUPPI jets, where Jet_jetId is not stored.
        # This implements a conservative tight-noise-ID style selection using
        # the available charged/neutral energy fractions and multiplicities.
        eta = abs(float(self.b_Jet_eta.At(i)))
        neHEF = float(self.b_Jet_neHEF.At(i))
        neEmEF = float(self.b_Jet_neEmEF.At(i))
        chHEF = float(self.b_Jet_chHEF.At(i))
        chEmEF = float(self.b_Jet_chEmEF.At(i))
        muEF = float(self.b_Jet_muEF.At(i))
        chMult = self._uchar_to_int(self.b_Jet_chMultiplicity.At(i))
        neMult = self._uchar_to_int(self.b_Jet_neMultiplicity.At(i))
        nConst = self._uchar_to_int(self.b_Jet_nConstituents.At(i))

        if eta <= 2.6:
            return (
                neHEF < 0.99 and
                neEmEF < 0.90 and
                nConst > 1 and
                muEF < 0.80 and
                chHEF > 0.01 and
                chMult > 0 and
                chEmEF < 0.80
            )
        if eta <= 2.7:
            return neHEF < 0.90 and neEmEF < 0.99 and muEF < 0.80
        if eta <= 3.0:
            return neEmEF < 0.99 and neMult > 2
        return neHEF < 0.90 and neEmEF < 0.90 and neMult > 10

    def isGoodJet(self, i):
        shifted_pt = self.b_Jet_pt.At(i)
        if self._is_jetvar_mode():
            shifted_pt = shifted_pt * (1.0 + self._jet_shift_sign() * self.b_Jet_Uncert.At(i))
        if not (shifted_pt > self.minJetPt_):
            return False
        if not (abs(self.b_Jet_eta.At(i)) < self.maxJetAbsEta_):
            return False
        official_jet_id = self._passes_official_ak4puppi_tight_jet_id(i)
        if official_jet_id is not None:
            if not official_jet_id:
                return False
        else:
            if not self._passes_tight_jet_id_from_composition(i):
                return False
        return True
    def analyze(self, event):
        if event._tree._ttreereaderversion > self._ttreereaderversion:
            self.initReaders(event._tree)

        self._reset_outputs()

        if self._is_jetvar_mode():
            return self._analyze_jetvar_only(event)

        self.out_MET_pt = float(self._event_value_first(event, ["PuppiMET_pt", "MET_T1_pt", "MET_pt", "PFMET_pt"]))
        self.out_MET_phi = float(self._event_value_first(event, ["PuppiMET_phi", "MET_T1_phi", "MET_phi", "PFMET_phi"]))

        muonIdxs = [i for i in range(self.b_Muon_pt.GetSize()) if self.isGoodMuon(i)]
        electronIdxs = [i for i in range(self.b_Electron_pt.GetSize()) if self.isGoodElectron(i)]
        self._sort_desc_by_pt(muonIdxs, self.b_Muon_pt)
        self._sort_desc_by_pt(electronIdxs, self.b_Electron_pt)

        self.out_nVetoMuon = len([i for i in range(self.b_Muon_pt.GetSize()) if self._passes_muon_veto(i)])
        self.out_nVetoElectron = len([i for i in range(self.b_Electron_pt.GetSize()) if self._passes_electron_veto(i)])

        npMuonIdxs = [i for i in range(self.b_Muon_pt.GetSize()) if self.isNPMuon(i)]
        npElectronIdxs = [i for i in range(self.b_Electron_pt.GetSize()) if self.isNPElectron(i)]
        self._sort_desc_by_pt(npMuonIdxs, self.b_Muon_pt)
        self._sort_desc_by_pt(npElectronIdxs, self.b_Electron_pt)

        tauIdxs = [i for i in range(self.b_Tau_pt.GetSize()) if self.isGoodTau(i, electronIdxs, muonIdxs)]
        npTauIdxs = [i for i in range(self.b_Tau_pt.GetSize()) if self.isNPTau(i, npElectronIdxs, npMuonIdxs)]
        self._sort_desc_by_pt(tauIdxs, self.b_Tau_pt)
        self._sort_desc_by_pt(npTauIdxs, self.b_Tau_pt)

        nGoodMuons = len(muonIdxs)
        nGoodElectrons = len(electronIdxs)
        nGoodTaus = len(tauIdxs)
        nNPMuons = len(npMuonIdxs)
        nNPElectrons = len(npElectronIdxs)
        nNPTaus = len(npTauIdxs)

        self.out_nGoodMuon = nGoodMuons
        self.out_nGoodElectron = nGoodElectrons
        self.out_nGoodTau = nGoodTaus
        self.out_nNPMuon = nNPMuons
        self.out_nNPElectron = nNPElectrons
        self.out_nNPTau = nNPTaus

        req_e, req_m, req_t = self._required_lepton_multiplicity(self.mode)
        if self.mode and self.mode.startswith("NPL"):
            if not (nNPElectrons == req_e and nNPMuons == req_m and nNPTaus == req_t):
                return False
        else:
            if not (nGoodElectrons == req_e and nGoodMuons == req_m and nGoodTaus == req_t):
                return False

        lepton1P4 = ROOT.TLorentzVector()
        lepton2P4 = ROOT.TLorentzVector()
        lepton3P4 = ROOT.TLorentzVector()
        lepton1jetP4 = ROOT.TLorentzVector()
        lepton2jetP4 = ROOT.TLorentzVector()
        lepton3jetP4 = ROOT.TLorentzVector()

        if nGoodTaus >= 1:
            i = tauIdxs[0]
            tauLeadP4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
            self.out_LeadingTau_p4 = self._copy_p4_from_lv(tauLeadP4)
        elif nNPTaus >= 1:
            i = npTauIdxs[0]
            tauLeadP4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
            self.out_LeadingTau_p4 = self._copy_p4_from_lv(tauLeadP4)

        # --------------------------------------------------------------
        # Central selections
        # --------------------------------------------------------------
        if self.mode == "MuElEl":
            if nGoodMuons >= 1 and nGoodElectrons >= 1:
                if nGoodMuons >= 1:
                    i = muonIdxs[0]
                    lepton1P4 = self._p4_from_arrays(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_phi, self.b_Muon_mass, i)
                    self.out_Lepton1_pdgId = -13 * int(self.b_Muon_charge.At(i))
                    self.out_LeadingMuon_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_PFIso = float(self.b_Muon_pfRelIso04_all.At(i))
                    lepton1jetP4 = self._p4_from_arrays(self.b_Muon_jetPt, self.b_Muon_jetEta, self.b_Muon_jetPhi, self.b_Muon_jetMass, i)
                    self.out_Lepton1_Jet_p4 = self._copy_p4_from_lv(lepton1jetP4)
                if nGoodElectrons >= 1:
                    i = electronIdxs[0]
                    lepton2P4 = self._p4_from_arrays(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_phi, self.b_Electron_mass, i)
                    self.out_Lepton2_pdgId = -11 * int(self.b_Electron_charge.At(i))
                    self.out_LeadingElectron_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_PFIso = float(self.b_Electron_pfRelIso03_all.At(i))
                    lepton2jetP4 = self._p4_from_arrays(self.b_Electron_jetPt, self.b_Electron_jetEta, self.b_Electron_jetPhi, self.b_Electron_jetMass, i)
                    self.out_Lepton2_Jet_p4 = self._copy_p4_from_lv(lepton2jetP4)
                if nGoodElectrons >= 2:
                    i = electronIdxs[1]
                    lepton3P4 = self._p4_from_arrays(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_phi, self.b_Electron_mass, i)
                    self.out_Lepton3_pdgId = -11 * int(self.b_Electron_charge.At(i))
                    self.out_Lepton3_p4 = self._copy_p4_from_lv(lepton3P4)
                    self.out_Lepton3_PFIso = float(self.b_Electron_pfRelIso03_all.At(i))
                    lepton3jetP4 = self._p4_from_arrays(self.b_Electron_jetPt, self.b_Electron_jetEta, self.b_Electron_jetPhi, self.b_Electron_jetMass, i)
                    self.out_Lepton3_Jet_p4 = self._copy_p4_from_lv(lepton3jetP4)

        elif self.mode == "ElMuMu":
            if nGoodElectrons >= 1 and nGoodMuons >= 1:
                if nGoodElectrons >= 1:
                    i = electronIdxs[0]
                    lepton1P4 = self._p4_from_arrays(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_phi, self.b_Electron_mass, i)
                    self.out_Lepton1_pdgId = -11 * int(self.b_Electron_charge.At(i))
                    self.out_LeadingElectron_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_PFIso = float(self.b_Electron_pfRelIso03_all.At(i))
                    lepton1jetP4 = self._p4_from_arrays(self.b_Electron_jetPt, self.b_Electron_jetEta, self.b_Electron_jetPhi, self.b_Electron_jetMass, i)
                    self.out_Lepton1_Jet_p4 = self._copy_p4_from_lv(lepton1jetP4)
                if nGoodMuons >= 1:
                    i = muonIdxs[0]
                    lepton2P4 = self._p4_from_arrays(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_phi, self.b_Muon_mass, i)
                    self.out_Lepton2_pdgId = -13 * int(self.b_Muon_charge.At(i))
                    self.out_LeadingMuon_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_PFIso = float(self.b_Muon_pfRelIso04_all.At(i))
                    lepton2jetP4 = self._p4_from_arrays(self.b_Muon_jetPt, self.b_Muon_jetEta, self.b_Muon_jetPhi, self.b_Muon_jetMass, i)
                    self.out_Lepton2_Jet_p4 = self._copy_p4_from_lv(lepton2jetP4)
                if nGoodMuons >= 2:
                    i = muonIdxs[1]
                    lepton3P4 = self._p4_from_arrays(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_phi, self.b_Muon_mass, i)
                    self.out_Lepton3_pdgId = -13 * int(self.b_Muon_charge.At(i))
                    # Intentional port of the original intermediate assignment typo.
                    self.out_Lepton3_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton3_PFIso = float(self.b_Muon_pfRelIso04_all.At(i))
                    lepton3jetP4 = self._p4_from_arrays(self.b_Muon_jetPt, self.b_Muon_jetEta, self.b_Muon_jetPhi, self.b_Muon_jetMass, i)
                    self.out_Lepton3_Jet_p4 = self._copy_p4_from_lv(lepton3jetP4)

        elif self.mode == "ElElEl":
            if nGoodElectrons >= 1:
                if nGoodElectrons >= 1:
                    i = electronIdxs[0]
                    lepton1P4 = self._p4_from_arrays(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_phi, self.b_Electron_mass, i)
                    self.out_Lepton1_pdgId = -11 * int(self.b_Electron_charge.At(i))
                    self.out_LeadingElectron_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_PFIso = float(self.b_Electron_pfRelIso03_all.At(i))
                    lepton1jetP4 = self._p4_from_arrays(self.b_Electron_jetPt, self.b_Electron_jetEta, self.b_Electron_jetPhi, self.b_Electron_jetMass, i)
                    self.out_Lepton1_Jet_p4 = self._copy_p4_from_lv(lepton1jetP4)
                if nGoodElectrons >= 2:
                    i = electronIdxs[1]
                    lepton2P4 = self._p4_from_arrays(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_phi, self.b_Electron_mass, i)
                    self.out_Lepton2_pdgId = -11 * int(self.b_Electron_charge.At(i))
                    self.out_Lepton2_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_PFIso = float(self.b_Electron_pfRelIso03_all.At(i))
                    lepton2jetP4 = self._p4_from_arrays(self.b_Electron_jetPt, self.b_Electron_jetEta, self.b_Electron_jetPhi, self.b_Electron_jetMass, i)
                    self.out_Lepton2_Jet_p4 = self._copy_p4_from_lv(lepton2jetP4)
                if nGoodElectrons >= 3:
                    i = electronIdxs[2]
                    lepton3P4 = self._p4_from_arrays(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_phi, self.b_Electron_mass, i)
                    self.out_Lepton3_pdgId = -11 * int(self.b_Electron_charge.At(i))
                    self.out_Lepton3_p4 = self._copy_p4_from_lv(lepton3P4)
                    self.out_Lepton3_PFIso = float(self.b_Electron_pfRelIso03_all.At(i))
                    lepton3jetP4 = self._p4_from_arrays(self.b_Electron_jetPt, self.b_Electron_jetEta, self.b_Electron_jetPhi, self.b_Electron_jetMass, i)
                    self.out_Lepton3_Jet_p4 = self._copy_p4_from_lv(lepton3jetP4)

        elif self.mode == "MuMuMu":
            if nGoodMuons >= 1:
                if nGoodMuons >= 1:
                    i = muonIdxs[0]
                    lepton1P4 = self._p4_from_arrays(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_phi, self.b_Muon_mass, i)
                    self.out_Lepton1_pdgId = -13 * int(self.b_Muon_charge.At(i))
                    self.out_LeadingMuon_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_PFIso = float(self.b_Muon_pfRelIso04_all.At(i))
                    lepton1jetP4 = self._p4_from_arrays(self.b_Muon_jetPt, self.b_Muon_jetEta, self.b_Muon_jetPhi, self.b_Muon_jetMass, i)
                    self.out_Lepton1_Jet_p4 = self._copy_p4_from_lv(lepton1jetP4)
                if nGoodMuons >= 2:
                    i = muonIdxs[1]
                    lepton2P4 = self._p4_from_arrays(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_phi, self.b_Muon_mass, i)
                    self.out_Lepton2_pdgId = -13 * int(self.b_Muon_charge.At(i))
                    self.out_Lepton2_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_PFIso = float(self.b_Muon_pfRelIso04_all.At(i))
                    lepton2jetP4 = self._p4_from_arrays(self.b_Muon_jetPt, self.b_Muon_jetEta, self.b_Muon_jetPhi, self.b_Muon_jetMass, i)
                    self.out_Lepton2_Jet_p4 = self._copy_p4_from_lv(lepton2jetP4)
                if nGoodMuons >= 3:
                    i = muonIdxs[2]
                    lepton3P4 = self._p4_from_arrays(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_phi, self.b_Muon_mass, i)
                    self.out_Lepton3_pdgId = -13 * int(self.b_Muon_charge.At(i))
                    self.out_Lepton3_p4 = self._copy_p4_from_lv(lepton3P4)
                    self.out_Lepton3_PFIso = float(self.b_Muon_pfRelIso04_all.At(i))
                    lepton3jetP4 = self._p4_from_arrays(self.b_Muon_jetPt, self.b_Muon_jetEta, self.b_Muon_jetPhi, self.b_Muon_jetMass, i)
                    self.out_Lepton3_Jet_p4 = self._copy_p4_from_lv(lepton3jetP4)

        elif self.mode == "TaTaTa":
            if nGoodTaus >= 1:
                if nGoodTaus >= 1:
                    i = tauIdxs[0]
                    lepton1P4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
                    self.out_Lepton1_pdgId = -15 * int(self.b_Tau_charge.At(i))
                    self.out_LeadingTau_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_PFIso = float(self.b_Tau_rawIso.At(i))
                    self.out_Lepton1_tauDecayMode = self._tau_decaymode(i)
                    self.out_Lepton1_tauGenPartFlav = self._tau_genpartflav(i)
                    self.out_Lepton1_tauProng = self._tau_prong_group_from_decaymode(self.out_Lepton1_tauDecayMode)
                    self.out_Lepton1_jetPtRatio = self._tau_jet_pt_ratio(i)
                    lepton1jetP4 = self._tau_matched_jet_p4(i)
                    self.out_Lepton1_Jet_p4 = self._copy_p4_from_lv(lepton1jetP4)
                    self.out_Lepton1_jetBTagDeepFlavor = self._tau_matched_jet_btag(i)
                if nGoodTaus >= 2:
                    i = tauIdxs[1]
                    lepton2P4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
                    self.out_Lepton2_pdgId = -15 * int(self.b_Tau_charge.At(i))
                    self.out_Lepton2_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_PFIso = float(self.b_Tau_rawIso.At(i))
                    self.out_Lepton2_tauDecayMode = self._tau_decaymode(i)
                    self.out_Lepton2_tauGenPartFlav = self._tau_genpartflav(i)
                    self.out_Lepton2_tauProng = self._tau_prong_group_from_decaymode(self.out_Lepton2_tauDecayMode)
                    self.out_Lepton2_jetPtRatio = self._tau_jet_pt_ratio(i)
                    lepton2jetP4 = self._tau_matched_jet_p4(i)
                    self.out_Lepton2_Jet_p4 = self._copy_p4_from_lv(lepton2jetP4)
                    self.out_Lepton2_jetBTagDeepFlavor = self._tau_matched_jet_btag(i)
                if nGoodTaus >= 3:
                    i = tauIdxs[2]
                    lepton3P4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
                    self.out_Lepton3_pdgId = -15 * int(self.b_Tau_charge.At(i))
                    self.out_Lepton3_p4 = self._copy_p4_from_lv(lepton3P4)
                    self.out_Lepton3_PFIso = float(self.b_Tau_rawIso.At(i))
                    self.out_Lepton3_tauDecayMode = self._tau_decaymode(i)
                    self.out_Lepton3_tauGenPartFlav = self._tau_genpartflav(i)
                    self.out_Lepton3_tauProng = self._tau_prong_group_from_decaymode(self.out_Lepton3_tauDecayMode)
                    self.out_Lepton3_jetPtRatio = self._tau_jet_pt_ratio(i)
                    lepton3jetP4 = self._tau_matched_jet_p4(i)
                    self.out_Lepton3_Jet_p4 = self._copy_p4_from_lv(lepton3jetP4)
                    self.out_Lepton3_jetBTagDeepFlavor = self._tau_matched_jet_btag(i)

        elif self.mode == "ElTaTa":
            if nGoodElectrons >= 1 and nGoodTaus >= 1:
                if nGoodElectrons >= 1:
                    i = electronIdxs[0]
                    lepton1P4 = self._p4_from_arrays(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_phi, self.b_Electron_mass, i)
                    self.out_Lepton1_pdgId = -11 * int(self.b_Electron_charge.At(i))
                    self.out_LeadingElectron_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_PFIso = float(self.b_Electron_pfRelIso03_all.At(i))
                    lepton1jetP4 = self._p4_from_arrays(self.b_Electron_jetPt, self.b_Electron_jetEta, self.b_Electron_jetPhi, self.b_Electron_jetMass, i)
                    self.out_Lepton1_Jet_p4 = self._copy_p4_from_lv(lepton1jetP4)
                if nGoodTaus >= 1:
                    i = tauIdxs[0]
                    lepton2P4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
                    self.out_Lepton2_pdgId = -15 * int(self.b_Tau_charge.At(i))
                    self.out_LeadingTau_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_PFIso = float(self.b_Tau_rawIso.At(i))
                    self.out_Lepton2_tauDecayMode = self._tau_decaymode(i)
                    self.out_Lepton2_tauGenPartFlav = self._tau_genpartflav(i)
                    self.out_Lepton2_tauProng = self._tau_prong_group_from_decaymode(self.out_Lepton2_tauDecayMode)
                    self.out_Lepton2_jetPtRatio = self._tau_jet_pt_ratio(i)
                    lepton2jetP4 = self._tau_matched_jet_p4(i)
                    self.out_Lepton2_Jet_p4 = self._copy_p4_from_lv(lepton2jetP4)
                    self.out_Lepton2_jetBTagDeepFlavor = self._tau_matched_jet_btag(i)
                if nGoodTaus >= 2:
                    i = tauIdxs[1]
                    lepton3P4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
                    self.out_Lepton3_pdgId = -15 * int(self.b_Tau_charge.At(i))
                    self.out_Lepton3_p4 = self._copy_p4_from_lv(lepton3P4)
                    self.out_Lepton3_PFIso = float(self.b_Tau_rawIso.At(i))
                    self.out_Lepton3_tauDecayMode = self._tau_decaymode(i)
                    self.out_Lepton3_tauGenPartFlav = self._tau_genpartflav(i)
                    self.out_Lepton3_tauProng = self._tau_prong_group_from_decaymode(self.out_Lepton3_tauDecayMode)
                    self.out_Lepton3_jetPtRatio = self._tau_jet_pt_ratio(i)
                    lepton3jetP4 = self._tau_matched_jet_p4(i)
                    self.out_Lepton3_Jet_p4 = self._copy_p4_from_lv(lepton3jetP4)
                    self.out_Lepton3_jetBTagDeepFlavor = self._tau_matched_jet_btag(i)

        elif self.mode == "MuTaTa":
            if nGoodMuons >= 1 and nGoodTaus >= 1:
                if nGoodMuons >= 1:
                    i = muonIdxs[0]
                    lepton1P4 = self._p4_from_arrays(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_phi, self.b_Muon_mass, i)
                    self.out_Lepton1_pdgId = -13 * int(self.b_Muon_charge.At(i))
                    self.out_LeadingMuon_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_PFIso = float(self.b_Muon_pfRelIso04_all.At(i))
                    lepton1jetP4 = self._p4_from_arrays(self.b_Muon_jetPt, self.b_Muon_jetEta, self.b_Muon_jetPhi, self.b_Muon_jetMass, i)
                    self.out_Lepton1_Jet_p4 = self._copy_p4_from_lv(lepton1jetP4)
                if nGoodTaus >= 1:
                    i = tauIdxs[0]
                    lepton2P4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
                    self.out_Lepton2_pdgId = -15 * int(self.b_Tau_charge.At(i))
                    self.out_LeadingTau_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_PFIso = float(self.b_Tau_rawIso.At(i))
                    self.out_Lepton2_tauDecayMode = self._tau_decaymode(i)
                    self.out_Lepton2_tauGenPartFlav = self._tau_genpartflav(i)
                    self.out_Lepton2_tauProng = self._tau_prong_group_from_decaymode(self.out_Lepton2_tauDecayMode)
                    self.out_Lepton2_jetPtRatio = self._tau_jet_pt_ratio(i)
                    lepton2jetP4 = self._tau_matched_jet_p4(i)
                    self.out_Lepton2_Jet_p4 = self._copy_p4_from_lv(lepton2jetP4)
                    self.out_Lepton2_jetBTagDeepFlavor = self._tau_matched_jet_btag(i)
                if nGoodTaus >= 2:
                    i = tauIdxs[1]
                    lepton3P4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
                    self.out_Lepton3_pdgId = -15 * int(self.b_Tau_charge.At(i))
                    self.out_Lepton3_p4 = self._copy_p4_from_lv(lepton3P4)
                    self.out_Lepton3_PFIso = float(self.b_Tau_rawIso.At(i))
                    self.out_Lepton3_tauDecayMode = self._tau_decaymode(i)
                    self.out_Lepton3_tauGenPartFlav = self._tau_genpartflav(i)
                    self.out_Lepton3_tauProng = self._tau_prong_group_from_decaymode(self.out_Lepton3_tauDecayMode)
                    self.out_Lepton3_jetPtRatio = self._tau_jet_pt_ratio(i)
                    lepton3jetP4 = self._tau_matched_jet_p4(i)
                    self.out_Lepton3_Jet_p4 = self._copy_p4_from_lv(lepton3jetP4)
                    self.out_Lepton3_jetBTagDeepFlavor = self._tau_matched_jet_btag(i)

        elif self.mode == "TaElEl":
            if nGoodTaus >= 1 and nGoodElectrons >= 1:
                if nGoodTaus >= 1:
                    i = tauIdxs[0]
                    lepton1P4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
                    self.out_Lepton1_pdgId = -15 * int(self.b_Tau_charge.At(i))
                    self.out_LeadingTau_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_PFIso = float(self.b_Tau_rawIso.At(i))
                    self.out_Lepton1_tauDecayMode = self._tau_decaymode(i)
                    self.out_Lepton1_tauGenPartFlav = self._tau_genpartflav(i)
                    self.out_Lepton1_tauProng = self._tau_prong_group_from_decaymode(self.out_Lepton1_tauDecayMode)
                    self.out_Lepton1_jetPtRatio = self._tau_jet_pt_ratio(i)
                    lepton1jetP4 = self._tau_matched_jet_p4(i)
                    self.out_Lepton1_Jet_p4 = self._copy_p4_from_lv(lepton1jetP4)
                    self.out_Lepton1_jetBTagDeepFlavor = self._tau_matched_jet_btag(i)
                if nGoodElectrons >= 1:
                    i = electronIdxs[0]
                    lepton2P4 = self._p4_from_arrays(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_phi, self.b_Electron_mass, i)
                    self.out_Lepton2_pdgId = -11 * int(self.b_Electron_charge.At(i))
                    self.out_LeadingElectron_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_PFIso = float(self.b_Electron_pfRelIso03_all.At(i))
                    lepton2jetP4 = self._p4_from_arrays(self.b_Electron_jetPt, self.b_Electron_jetEta, self.b_Electron_jetPhi, self.b_Electron_jetMass, i)
                    self.out_Lepton2_Jet_p4 = self._copy_p4_from_lv(lepton2jetP4)
                if nGoodElectrons >= 2:
                    i = electronIdxs[1]
                    lepton3P4 = self._p4_from_arrays(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_phi, self.b_Electron_mass, i)
                    self.out_Lepton3_pdgId = -11 * int(self.b_Electron_charge.At(i))
                    self.out_Lepton3_p4 = self._copy_p4_from_lv(lepton3P4)
                    self.out_Lepton3_PFIso = float(self.b_Electron_pfRelIso03_all.At(i))
                    lepton3jetP4 = self._p4_from_arrays(self.b_Electron_jetPt, self.b_Electron_jetEta, self.b_Electron_jetPhi, self.b_Electron_jetMass, i)
                    self.out_Lepton3_Jet_p4 = self._copy_p4_from_lv(lepton3jetP4)

        elif self.mode == "TaMuMu":
            if nGoodTaus >= 1 and nGoodMuons >= 1:
                if nGoodTaus >= 1:
                    i = tauIdxs[0]
                    lepton1P4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
                    self.out_Lepton1_pdgId = -15 * int(self.b_Tau_charge.At(i))
                    self.out_LeadingTau_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_PFIso = float(self.b_Tau_rawIso.At(i))
                    self.out_Lepton1_tauDecayMode = self._tau_decaymode(i)
                    self.out_Lepton1_tauGenPartFlav = self._tau_genpartflav(i)
                    self.out_Lepton1_tauProng = self._tau_prong_group_from_decaymode(self.out_Lepton1_tauDecayMode)
                    self.out_Lepton1_jetPtRatio = self._tau_jet_pt_ratio(i)
                    lepton1jetP4 = self._tau_matched_jet_p4(i)
                    self.out_Lepton1_Jet_p4 = self._copy_p4_from_lv(lepton1jetP4)
                    self.out_Lepton1_jetBTagDeepFlavor = self._tau_matched_jet_btag(i)
                if nGoodMuons >= 1:
                    i = muonIdxs[0]
                    lepton2P4 = self._p4_from_arrays(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_phi, self.b_Muon_mass, i)
                    self.out_Lepton2_pdgId = -13 * int(self.b_Muon_charge.At(i))
                    self.out_LeadingMuon_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_PFIso = float(self.b_Muon_pfRelIso04_all.At(i))
                    lepton2jetP4 = self._p4_from_arrays(self.b_Muon_jetPt, self.b_Muon_jetEta, self.b_Muon_jetPhi, self.b_Muon_jetMass, i)
                    self.out_Lepton2_Jet_p4 = self._copy_p4_from_lv(lepton2jetP4)
                if nGoodMuons >= 2:
                    i = muonIdxs[1]
                    lepton3P4 = self._p4_from_arrays(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_phi, self.b_Muon_mass, i)
                    self.out_Lepton3_pdgId = -13 * int(self.b_Muon_charge.At(i))
                    self.out_Lepton3_p4 = self._copy_p4_from_lv(lepton3P4)
                    self.out_Lepton3_PFIso = float(self.b_Muon_pfRelIso04_all.At(i))
                    lepton3jetP4 = self._p4_from_arrays(self.b_Muon_jetPt, self.b_Muon_jetEta, self.b_Muon_jetPhi, self.b_Muon_jetMass, i)
                    self.out_Lepton3_Jet_p4 = self._copy_p4_from_lv(lepton3jetP4)

        # --------------------------------------------------------------
        # NPL selections
        # --------------------------------------------------------------
        elif self.mode == "NPLMuElEl":
            if nNPMuons >= 1 and nNPElectrons >= 1:
                if nNPMuons >= 1:
                    i = npMuonIdxs[0]
                    lepton1P4 = self._p4_from_arrays(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_phi, self.b_Muon_mass, i)
                    self.out_Lepton1_pdgId = -13 * int(self.b_Muon_charge.At(i))
                    self.out_LeadingMuon_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_PFIso = float(self.b_Muon_pfRelIso04_all.At(i))
                    self.out_Lepton1_jetPtRatio = float(self.b_Muon_jetPtRatio.At(i))
                    lepton1jetP4 = self._p4_from_arrays(self.b_Muon_jetPt, self.b_Muon_jetEta, self.b_Muon_jetPhi, self.b_Muon_jetMass, i)
                    self.out_Lepton1_Jet_p4 = self._copy_p4_from_lv(lepton1jetP4)
                    self.out_Lepton1_mvaTOP = float(self.b_Muon_mvaTOP.At(i))
                    if self.isGoodMuon(i):
                        self.out_Lepton1_isTight = True
                    if self.isOnlyFakeMuon(i):
                        self.out_Lepton1_isFake = True
                if nNPElectrons >= 1:
                    i = npElectronIdxs[0]
                    lepton2P4 = self._p4_from_arrays(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_phi, self.b_Electron_mass, i)
                    self.out_Lepton2_pdgId = -11 * int(self.b_Electron_charge.At(i))
                    self.out_LeadingElectron_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_PFIso = float(self.b_Electron_pfRelIso03_all.At(i))
                    self.out_Lepton2_jetPtRatio = float(self.b_Electron_jetPtRatio.At(i))
                    lepton2jetP4 = self._p4_from_arrays(self.b_Electron_jetPt, self.b_Electron_jetEta, self.b_Electron_jetPhi, self.b_Electron_jetMass, i)
                    self.out_Lepton2_Jet_p4 = self._copy_p4_from_lv(lepton2jetP4)
                    self.out_Lepton2_mvaTOP = float(self.b_Electron_mvaTOP.At(i))
                    if self.isGoodElectron(i):
                        self.out_Lepton2_isTight = True
                    if self.isOnlyFakeElectron(i):
                        self.out_Lepton2_isFake = True
                if nNPElectrons >= 2:
                    i = npElectronIdxs[1]
                    lepton3P4 = self._p4_from_arrays(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_phi, self.b_Electron_mass, i)
                    self.out_Lepton3_pdgId = -11 * int(self.b_Electron_charge.At(i))
                    self.out_Lepton3_p4 = self._copy_p4_from_lv(lepton3P4)
                    self.out_Lepton3_PFIso = float(self.b_Electron_pfRelIso03_all.At(i))
                    self.out_Lepton3_jetPtRatio = float(self.b_Electron_jetPtRatio.At(i))
                    lepton3jetP4 = self._p4_from_arrays(self.b_Electron_jetPt, self.b_Electron_jetEta, self.b_Electron_jetPhi, self.b_Electron_jetMass, i)
                    self.out_Lepton3_Jet_p4 = self._copy_p4_from_lv(lepton3jetP4)
                    # Intentional port of the original typo in FCNCTriLeptonCppWorker.cc.
                    self.out_Lepton2_mvaTOP = float(self.b_Electron_mvaTOP.At(i))
                    if self.isGoodElectron(i):
                        self.out_Lepton3_isTight = True
                    if self.isOnlyFakeElectron(i):
                        self.out_Lepton3_isFake = True

        elif self.mode == "NPLElMuMu":
            if nNPElectrons >= 1 and nNPMuons >= 1:
                if nNPElectrons >= 1:
                    i = npElectronIdxs[0]
                    lepton1P4 = self._p4_from_arrays(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_phi, self.b_Electron_mass, i)
                    self.out_Lepton1_pdgId = -11 * int(self.b_Electron_charge.At(i))
                    self.out_LeadingElectron_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_PFIso = float(self.b_Electron_pfRelIso03_all.At(i))
                    self.out_Lepton1_jetPtRatio = float(self.b_Electron_jetPtRatio.At(i))
                    lepton1jetP4 = self._p4_from_arrays(self.b_Electron_jetPt, self.b_Electron_jetEta, self.b_Electron_jetPhi, self.b_Electron_jetMass, i)
                    self.out_Lepton1_Jet_p4 = self._copy_p4_from_lv(lepton1jetP4)
                    self.out_Lepton1_mvaTOP = float(self.b_Electron_mvaTOP.At(i))
                    if self.isGoodElectron(i):
                        self.out_Lepton1_isTight = True
                    if self.isOnlyFakeElectron(i):
                        self.out_Lepton1_isFake = True
                if nNPMuons >= 1:
                    i = npMuonIdxs[0]
                    lepton2P4 = self._p4_from_arrays(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_phi, self.b_Muon_mass, i)
                    self.out_Lepton2_pdgId = -13 * int(self.b_Muon_charge.At(i))
                    self.out_LeadingMuon_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_PFIso = float(self.b_Muon_pfRelIso04_all.At(i))
                    self.out_Lepton2_jetPtRatio = float(self.b_Muon_jetPtRatio.At(i))
                    lepton2jetP4 = self._p4_from_arrays(self.b_Muon_jetPt, self.b_Muon_jetEta, self.b_Muon_jetPhi, self.b_Muon_jetMass, i)
                    self.out_Lepton2_Jet_p4 = self._copy_p4_from_lv(lepton2jetP4)
                    self.out_Lepton2_mvaTOP = float(self.b_Muon_mvaTOP.At(i))
                    if self.isGoodMuon(i):
                        self.out_Lepton2_isTight = True
                    if self.isOnlyFakeMuon(i):
                        self.out_Lepton2_isFake = True
                if nNPMuons >= 2:
                    i = npMuonIdxs[1]
                    lepton3P4 = self._p4_from_arrays(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_phi, self.b_Muon_mass, i)
                    self.out_Lepton3_pdgId = -13 * int(self.b_Muon_charge.At(i))
                    self.out_Lepton3_p4 = self._copy_p4_from_lv(lepton3P4)
                    self.out_Lepton3_PFIso = float(self.b_Muon_pfRelIso04_all.At(i))
                    self.out_Lepton3_jetPtRatio = float(self.b_Muon_jetPtRatio.At(i))
                    lepton3jetP4 = self._p4_from_arrays(self.b_Muon_jetPt, self.b_Muon_jetEta, self.b_Muon_jetPhi, self.b_Muon_jetMass, i)
                    self.out_Lepton3_Jet_p4 = self._copy_p4_from_lv(lepton3jetP4)
                    self.out_Lepton3_mvaTOP = float(self.b_Muon_mvaTOP.At(i))
                    if self.isGoodMuon(i):
                        self.out_Lepton3_isTight = True
                    if self.isOnlyFakeMuon(i):
                        self.out_Lepton3_isFake = True

        elif self.mode == "NPLElElEl":
            if nNPElectrons >= 1:
                if nNPElectrons >= 1:
                    i = npElectronIdxs[0]
                    lepton1P4 = self._p4_from_arrays(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_phi, self.b_Electron_mass, i)
                    self.out_Lepton1_pdgId = -11 * int(self.b_Electron_charge.At(i))
                    self.out_LeadingElectron_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_PFIso = float(self.b_Electron_pfRelIso03_all.At(i))
                    self.out_Lepton1_jetPtRatio = float(self.b_Electron_jetPtRatio.At(i))
                    lepton1jetP4 = self._p4_from_arrays(self.b_Electron_jetPt, self.b_Electron_jetEta, self.b_Electron_jetPhi, self.b_Electron_jetMass, i)
                    self.out_Lepton1_Jet_p4 = self._copy_p4_from_lv(lepton1jetP4)
                    self.out_Lepton1_mvaTOP = float(self.b_Electron_mvaTOP.At(i))
                    if self.isGoodElectron(i):
                        self.out_Lepton1_isTight = True
                    if self.isOnlyFakeElectron(i):
                        self.out_Lepton1_isFake = True
                if nNPElectrons >= 2:
                    i = npElectronIdxs[1]
                    lepton2P4 = self._p4_from_arrays(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_phi, self.b_Electron_mass, i)
                    self.out_Lepton2_pdgId = -11 * int(self.b_Electron_charge.At(i))
                    self.out_Lepton2_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_PFIso = float(self.b_Electron_pfRelIso03_all.At(i))
                    self.out_Lepton2_jetPtRatio = float(self.b_Electron_jetPtRatio.At(i))
                    lepton2jetP4 = self._p4_from_arrays(self.b_Electron_jetPt, self.b_Electron_jetEta, self.b_Electron_jetPhi, self.b_Electron_jetMass, i)
                    self.out_Lepton2_Jet_p4 = self._copy_p4_from_lv(lepton2jetP4)
                    self.out_Lepton2_mvaTOP = float(self.b_Electron_mvaTOP.At(i))
                    if self.isGoodElectron(i):
                        self.out_Lepton2_isTight = True
                    if self.isOnlyFakeElectron(i):
                        self.out_Lepton2_isFake = True
                if nNPElectrons >= 3:
                    i = npElectronIdxs[2]
                    lepton3P4 = self._p4_from_arrays(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_phi, self.b_Electron_mass, i)
                    self.out_Lepton3_pdgId = -11 * int(self.b_Electron_charge.At(i))
                    self.out_Lepton3_p4 = self._copy_p4_from_lv(lepton3P4)
                    self.out_Lepton3_PFIso = float(self.b_Electron_pfRelIso03_all.At(i))
                    self.out_Lepton3_jetPtRatio = float(self.b_Electron_jetPtRatio.At(i))
                    lepton3jetP4 = self._p4_from_arrays(self.b_Electron_jetPt, self.b_Electron_jetEta, self.b_Electron_jetPhi, self.b_Electron_jetMass, i)
                    self.out_Lepton3_Jet_p4 = self._copy_p4_from_lv(lepton3jetP4)
                    self.out_Lepton3_mvaTOP = float(self.b_Electron_mvaTOP.At(i))
                    if self.isGoodElectron(i):
                        self.out_Lepton3_isTight = True
                    if self.isOnlyFakeElectron(i):
                        self.out_Lepton3_isFake = True

        elif self.mode == "NPLMuMuMu":
            if nNPMuons >= 1:
                if nNPMuons >= 1:
                    i = npMuonIdxs[0]
                    lepton1P4 = self._p4_from_arrays(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_phi, self.b_Muon_mass, i)
                    self.out_Lepton1_pdgId = -13 * int(self.b_Muon_charge.At(i))
                    self.out_LeadingMuon_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_PFIso = float(self.b_Muon_pfRelIso04_all.At(i))
                    self.out_Lepton1_jetPtRatio = float(self.b_Muon_jetPtRatio.At(i))
                    lepton1jetP4 = self._p4_from_arrays(self.b_Muon_jetPt, self.b_Muon_jetEta, self.b_Muon_jetPhi, self.b_Muon_jetMass, i)
                    self.out_Lepton1_Jet_p4 = self._copy_p4_from_lv(lepton1jetP4)
                    self.out_Lepton1_mvaTOP = float(self.b_Muon_mvaTOP.At(i))
                    if self.isGoodMuon(i):
                        self.out_Lepton1_isTight = True
                    if self.isOnlyFakeMuon(i):
                        self.out_Lepton1_isFake = True
                if nNPMuons >= 2:
                    i = npMuonIdxs[1]
                    lepton2P4 = self._p4_from_arrays(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_phi, self.b_Muon_mass, i)
                    self.out_Lepton2_pdgId = -13 * int(self.b_Muon_charge.At(i))
                    self.out_Lepton2_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_PFIso = float(self.b_Muon_pfRelIso04_all.At(i))
                    self.out_Lepton2_jetPtRatio = float(self.b_Muon_jetPtRatio.At(i))
                    lepton2jetP4 = self._p4_from_arrays(self.b_Muon_jetPt, self.b_Muon_jetEta, self.b_Muon_jetPhi, self.b_Muon_jetMass, i)
                    self.out_Lepton2_Jet_p4 = self._copy_p4_from_lv(lepton2jetP4)
                    self.out_Lepton2_mvaTOP = float(self.b_Muon_mvaTOP.At(i))
                    if self.isGoodMuon(i):
                        self.out_Lepton2_isTight = True
                    if self.isOnlyFakeMuon(i):
                        self.out_Lepton2_isFake = True
                if nNPMuons >= 3:
                    i = npMuonIdxs[2]
                    lepton3P4 = self._p4_from_arrays(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_phi, self.b_Muon_mass, i)
                    self.out_Lepton3_pdgId = -13 * int(self.b_Muon_charge.At(i))
                    self.out_Lepton3_p4 = self._copy_p4_from_lv(lepton3P4)
                    self.out_Lepton3_PFIso = float(self.b_Muon_pfRelIso04_all.At(i))
                    self.out_Lepton3_jetPtRatio = float(self.b_Muon_jetPtRatio.At(i))
                    lepton3jetP4 = self._p4_from_arrays(self.b_Muon_jetPt, self.b_Muon_jetEta, self.b_Muon_jetPhi, self.b_Muon_jetMass, i)
                    self.out_Lepton3_Jet_p4 = self._copy_p4_from_lv(lepton3jetP4)
                    self.out_Lepton3_mvaTOP = float(self.b_Muon_mvaTOP.At(i))
                    if self.isGoodMuon(i):
                        self.out_Lepton3_isTight = True
                    if self.isOnlyFakeMuon(i):
                        self.out_Lepton3_isFake = True

        elif self.mode == "NPLTaTaTa":
            if nNPTaus >= 1:
                if nNPTaus >= 1:
                    i = npTauIdxs[0]
                    lepton1P4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
                    self.out_Lepton1_pdgId = -15 * int(self.b_Tau_charge.At(i))
                    self.out_LeadingTau_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_PFIso = float(self.b_Tau_rawIso.At(i))
                    self.out_Lepton1_tauDecayMode = self._tau_decaymode(i)
                    self.out_Lepton1_tauGenPartFlav = self._tau_genpartflav(i)
                    self.out_Lepton1_tauProng = self._tau_prong_group_from_decaymode(self.out_Lepton1_tauDecayMode)
                    self.out_Lepton1_jetPtRatio = self._tau_jet_pt_ratio(i)
                    lepton1jetP4 = self._tau_matched_jet_p4(i)
                    self.out_Lepton1_Jet_p4 = self._copy_p4_from_lv(lepton1jetP4)
                    self.out_Lepton1_jetBTagDeepFlavor = self._tau_matched_jet_btag(i)
                    if self.isGoodTau(i, npElectronIdxs, npMuonIdxs):
                        self.out_Lepton1_isTight = True
                    if self.isOnlyFakeTau(i, npElectronIdxs, npMuonIdxs):
                        self.out_Lepton1_isFake = True
                if nNPTaus >= 2:
                    i = npTauIdxs[1]
                    lepton2P4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
                    self.out_Lepton2_pdgId = -15 * int(self.b_Tau_charge.At(i))
                    self.out_Lepton2_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_PFIso = float(self.b_Tau_rawIso.At(i))
                    self.out_Lepton2_tauDecayMode = self._tau_decaymode(i)
                    self.out_Lepton2_tauGenPartFlav = self._tau_genpartflav(i)
                    self.out_Lepton2_tauProng = self._tau_prong_group_from_decaymode(self.out_Lepton2_tauDecayMode)
                    self.out_Lepton2_jetPtRatio = self._tau_jet_pt_ratio(i)
                    lepton2jetP4 = self._tau_matched_jet_p4(i)
                    self.out_Lepton2_Jet_p4 = self._copy_p4_from_lv(lepton2jetP4)
                    self.out_Lepton2_jetBTagDeepFlavor = self._tau_matched_jet_btag(i)
                    if self.isGoodTau(i, npElectronIdxs, npMuonIdxs):
                        self.out_Lepton2_isTight = True
                    if self.isOnlyFakeTau(i, npElectronIdxs, npMuonIdxs):
                        self.out_Lepton2_isFake = True
                if nNPTaus >= 3:
                    i = npTauIdxs[2]
                    lepton3P4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
                    self.out_Lepton3_pdgId = -15 * int(self.b_Tau_charge.At(i))
                    self.out_Lepton3_p4 = self._copy_p4_from_lv(lepton3P4)
                    self.out_Lepton3_PFIso = float(self.b_Tau_rawIso.At(i))
                    self.out_Lepton3_tauDecayMode = self._tau_decaymode(i)
                    self.out_Lepton3_tauGenPartFlav = self._tau_genpartflav(i)
                    self.out_Lepton3_tauProng = self._tau_prong_group_from_decaymode(self.out_Lepton3_tauDecayMode)
                    self.out_Lepton3_jetPtRatio = self._tau_jet_pt_ratio(i)
                    lepton3jetP4 = self._tau_matched_jet_p4(i)
                    self.out_Lepton3_Jet_p4 = self._copy_p4_from_lv(lepton3jetP4)
                    self.out_Lepton3_jetBTagDeepFlavor = self._tau_matched_jet_btag(i)
                    if self.isGoodTau(i, npElectronIdxs, npMuonIdxs):
                        self.out_Lepton3_isTight = True
                    if self.isOnlyFakeTau(i, npElectronIdxs, npMuonIdxs):
                        self.out_Lepton3_isFake = True

        elif self.mode == "NPLElTaTa":
            if nNPElectrons >= 1 and nNPTaus >= 1:
                if nNPElectrons >= 1:
                    i = npElectronIdxs[0]
                    lepton1P4 = self._p4_from_arrays(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_phi, self.b_Electron_mass, i)
                    self.out_Lepton1_pdgId = -11 * int(self.b_Electron_charge.At(i))
                    self.out_LeadingElectron_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_PFIso = float(self.b_Electron_pfRelIso03_all.At(i))
                    self.out_Lepton1_jetPtRatio = float(self.b_Electron_jetPtRatio.At(i))
                    lepton1jetP4 = self._p4_from_arrays(self.b_Electron_jetPt, self.b_Electron_jetEta, self.b_Electron_jetPhi, self.b_Electron_jetMass, i)
                    self.out_Lepton1_Jet_p4 = self._copy_p4_from_lv(lepton1jetP4)
                    self.out_Lepton1_mvaTOP = float(self.b_Electron_mvaTOP.At(i))
                    if self.isGoodElectron(i):
                        self.out_Lepton1_isTight = True
                    if self.isOnlyFakeElectron(i):
                        self.out_Lepton1_isFake = True
                if nNPTaus >= 1:
                    i = npTauIdxs[0]
                    lepton2P4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
                    self.out_Lepton2_pdgId = -15 * int(self.b_Tau_charge.At(i))
                    self.out_Lepton2_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_PFIso = float(self.b_Tau_rawIso.At(i))
                    self.out_Lepton2_tauDecayMode = self._tau_decaymode(i)
                    self.out_Lepton2_tauGenPartFlav = self._tau_genpartflav(i)
                    self.out_Lepton2_tauProng = self._tau_prong_group_from_decaymode(self.out_Lepton2_tauDecayMode)
                    self.out_Lepton2_jetPtRatio = self._tau_jet_pt_ratio(i)
                    lepton2jetP4 = self._tau_matched_jet_p4(i)
                    self.out_Lepton2_Jet_p4 = self._copy_p4_from_lv(lepton2jetP4)
                    self.out_Lepton2_jetBTagDeepFlavor = self._tau_matched_jet_btag(i)
                    if self.isGoodTau(i, npElectronIdxs, npMuonIdxs):
                        self.out_Lepton2_isTight = True
                    if self.isOnlyFakeTau(i, npElectronIdxs, npMuonIdxs):
                        self.out_Lepton2_isFake = True
                if nNPTaus >= 2:
                    i = npTauIdxs[1]
                    lepton3P4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
                    self.out_Lepton3_pdgId = -15 * int(self.b_Tau_charge.At(i))
                    self.out_Lepton3_p4 = self._copy_p4_from_lv(lepton3P4)
                    self.out_Lepton3_PFIso = float(self.b_Tau_rawIso.At(i))
                    self.out_Lepton3_tauDecayMode = self._tau_decaymode(i)
                    self.out_Lepton3_tauGenPartFlav = self._tau_genpartflav(i)
                    self.out_Lepton3_tauProng = self._tau_prong_group_from_decaymode(self.out_Lepton3_tauDecayMode)
                    self.out_Lepton3_jetPtRatio = self._tau_jet_pt_ratio(i)
                    lepton3jetP4 = self._tau_matched_jet_p4(i)
                    self.out_Lepton3_Jet_p4 = self._copy_p4_from_lv(lepton3jetP4)
                    self.out_Lepton3_jetBTagDeepFlavor = self._tau_matched_jet_btag(i)
                    if self.isGoodTau(i, npElectronIdxs, npMuonIdxs):
                        self.out_Lepton3_isTight = True
                    if self.isOnlyFakeTau(i, npElectronIdxs, npMuonIdxs):
                        self.out_Lepton3_isFake = True

        elif self.mode == "NPLMuTaTa":
            if nNPMuons >= 1 and nNPTaus >= 1:
                if nNPMuons >= 1:
                    i = npMuonIdxs[0]
                    lepton1P4 = self._p4_from_arrays(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_phi, self.b_Muon_mass, i)
                    self.out_Lepton1_pdgId = -13 * int(self.b_Muon_charge.At(i))
                    self.out_LeadingMuon_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_PFIso = float(self.b_Muon_pfRelIso04_all.At(i))
                    self.out_Lepton1_jetPtRatio = float(self.b_Muon_jetPtRatio.At(i))
                    lepton1jetP4 = self._p4_from_arrays(self.b_Muon_jetPt, self.b_Muon_jetEta, self.b_Muon_jetPhi, self.b_Muon_jetMass, i)
                    self.out_Lepton1_Jet_p4 = self._copy_p4_from_lv(lepton1jetP4)
                    self.out_Lepton1_mvaTOP = float(self.b_Muon_mvaTOP.At(i))
                    if self.isGoodMuon(i):
                        self.out_Lepton1_isTight = True
                    if self.isOnlyFakeMuon(i):
                        self.out_Lepton1_isFake = True
                if nNPTaus >= 1:
                    i = npTauIdxs[0]
                    lepton2P4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
                    self.out_Lepton2_pdgId = -15 * int(self.b_Tau_charge.At(i))
                    self.out_Lepton2_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_PFIso = float(self.b_Tau_rawIso.At(i))
                    self.out_Lepton2_tauDecayMode = self._tau_decaymode(i)
                    self.out_Lepton2_tauGenPartFlav = self._tau_genpartflav(i)
                    self.out_Lepton2_tauProng = self._tau_prong_group_from_decaymode(self.out_Lepton2_tauDecayMode)
                    self.out_Lepton2_jetPtRatio = self._tau_jet_pt_ratio(i)
                    lepton2jetP4 = self._tau_matched_jet_p4(i)
                    self.out_Lepton2_Jet_p4 = self._copy_p4_from_lv(lepton2jetP4)
                    self.out_Lepton2_jetBTagDeepFlavor = self._tau_matched_jet_btag(i)
                    if self.isGoodTau(i, npElectronIdxs, npMuonIdxs):
                        self.out_Lepton2_isTight = True
                    if self.isOnlyFakeTau(i, npElectronIdxs, npMuonIdxs):
                        self.out_Lepton2_isFake = True
                if nNPTaus >= 2:
                    i = npTauIdxs[1]
                    lepton3P4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
                    self.out_Lepton3_pdgId = -15 * int(self.b_Tau_charge.At(i))
                    self.out_Lepton3_p4 = self._copy_p4_from_lv(lepton3P4)
                    self.out_Lepton3_PFIso = float(self.b_Tau_rawIso.At(i))
                    self.out_Lepton3_tauDecayMode = self._tau_decaymode(i)
                    self.out_Lepton3_tauGenPartFlav = self._tau_genpartflav(i)
                    self.out_Lepton3_tauProng = self._tau_prong_group_from_decaymode(self.out_Lepton3_tauDecayMode)
                    self.out_Lepton3_jetPtRatio = self._tau_jet_pt_ratio(i)
                    lepton3jetP4 = self._tau_matched_jet_p4(i)
                    self.out_Lepton3_Jet_p4 = self._copy_p4_from_lv(lepton3jetP4)
                    self.out_Lepton3_jetBTagDeepFlavor = self._tau_matched_jet_btag(i)
                    if self.isGoodTau(i, npElectronIdxs, npMuonIdxs):
                        self.out_Lepton3_isTight = True
                    if self.isOnlyFakeTau(i, npElectronIdxs, npMuonIdxs):
                        self.out_Lepton3_isFake = True

        elif self.mode == "NPLTaElEl":
            if nNPTaus >= 1 and nNPElectrons >= 1:
                if nNPTaus >= 1:
                    i = npTauIdxs[0]
                    lepton1P4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
                    self.out_Lepton1_pdgId = -15 * int(self.b_Tau_charge.At(i))
                    self.out_LeadingTau_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_PFIso = float(self.b_Tau_rawIso.At(i))
                    self.out_Lepton1_tauDecayMode = self._tau_decaymode(i)
                    self.out_Lepton1_tauGenPartFlav = self._tau_genpartflav(i)
                    self.out_Lepton1_tauProng = self._tau_prong_group_from_decaymode(self.out_Lepton1_tauDecayMode)
                    self.out_Lepton1_jetPtRatio = self._tau_jet_pt_ratio(i)
                    lepton1jetP4 = self._tau_matched_jet_p4(i)
                    self.out_Lepton1_Jet_p4 = self._copy_p4_from_lv(lepton1jetP4)
                    self.out_Lepton1_jetBTagDeepFlavor = self._tau_matched_jet_btag(i)
                    if self.isGoodTau(i, npElectronIdxs, npMuonIdxs):
                        self.out_Lepton1_isTight = True
                    if self.isOnlyFakeTau(i, npElectronIdxs, npMuonIdxs):
                        self.out_Lepton1_isFake = True
                if nNPElectrons >= 1:
                    i = npElectronIdxs[0]
                    lepton2P4 = self._p4_from_arrays(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_phi, self.b_Electron_mass, i)
                    self.out_Lepton2_pdgId = -11 * int(self.b_Electron_charge.At(i))
                    self.out_LeadingElectron_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_PFIso = float(self.b_Electron_pfRelIso03_all.At(i))
                    self.out_Lepton2_jetPtRatio = float(self.b_Electron_jetPtRatio.At(i))
                    lepton2jetP4 = self._p4_from_arrays(self.b_Electron_jetPt, self.b_Electron_jetEta, self.b_Electron_jetPhi, self.b_Electron_jetMass, i)
                    self.out_Lepton2_Jet_p4 = self._copy_p4_from_lv(lepton2jetP4)
                    self.out_Lepton2_mvaTOP = float(self.b_Electron_mvaTOP.At(i))
                    if self.isGoodElectron(i):
                        self.out_Lepton2_isTight = True
                    if self.isOnlyFakeElectron(i):
                        self.out_Lepton2_isFake = True
                if nNPElectrons >= 2:
                    i = npElectronIdxs[1]
                    lepton3P4 = self._p4_from_arrays(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_phi, self.b_Electron_mass, i)
                    self.out_Lepton3_pdgId = -11 * int(self.b_Electron_charge.At(i))
                    self.out_Lepton3_p4 = self._copy_p4_from_lv(lepton3P4)
                    self.out_Lepton3_PFIso = float(self.b_Electron_pfRelIso03_all.At(i))
                    self.out_Lepton3_jetPtRatio = float(self.b_Electron_jetPtRatio.At(i))
                    lepton3jetP4 = self._p4_from_arrays(self.b_Electron_jetPt, self.b_Electron_jetEta, self.b_Electron_jetPhi, self.b_Electron_jetMass, i)
                    self.out_Lepton3_Jet_p4 = self._copy_p4_from_lv(lepton3jetP4)
                    self.out_Lepton3_mvaTOP = float(self.b_Electron_mvaTOP.At(i))
                    if self.isGoodElectron(i):
                        self.out_Lepton3_isTight = True
                    if self.isOnlyFakeElectron(i):
                        self.out_Lepton3_isFake = True

        elif self.mode == "NPLTaMuMu":
            if nNPTaus >= 1 and nNPMuons >= 1:
                if nNPTaus >= 1:
                    i = npTauIdxs[0]
                    lepton1P4 = self._p4_from_arrays(self.b_Tau_pt, self.b_Tau_eta, self.b_Tau_phi, self.b_Tau_mass, i)
                    self.out_Lepton1_pdgId = -15 * int(self.b_Tau_charge.At(i))
                    self.out_LeadingTau_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_p4 = self._copy_p4_from_lv(lepton1P4)
                    self.out_Lepton1_PFIso = float(self.b_Tau_rawIso.At(i))
                    self.out_Lepton1_tauDecayMode = self._tau_decaymode(i)
                    self.out_Lepton1_tauGenPartFlav = self._tau_genpartflav(i)
                    self.out_Lepton1_tauProng = self._tau_prong_group_from_decaymode(self.out_Lepton1_tauDecayMode)
                    self.out_Lepton1_jetPtRatio = self._tau_jet_pt_ratio(i)
                    lepton1jetP4 = self._tau_matched_jet_p4(i)
                    self.out_Lepton1_Jet_p4 = self._copy_p4_from_lv(lepton1jetP4)
                    self.out_Lepton1_jetBTagDeepFlavor = self._tau_matched_jet_btag(i)
                    if self.isGoodTau(i, npElectronIdxs, npMuonIdxs):
                        self.out_Lepton1_isTight = True
                    if self.isOnlyFakeTau(i, npElectronIdxs, npMuonIdxs):
                        self.out_Lepton1_isFake = True
                if nNPMuons >= 1:
                    i = npMuonIdxs[0]
                    lepton2P4 = self._p4_from_arrays(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_phi, self.b_Muon_mass, i)
                    self.out_Lepton2_pdgId = -13 * int(self.b_Muon_charge.At(i))
                    self.out_LeadingMuon_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_p4 = self._copy_p4_from_lv(lepton2P4)
                    self.out_Lepton2_PFIso = float(self.b_Muon_pfRelIso04_all.At(i))
                    self.out_Lepton2_jetPtRatio = float(self.b_Muon_jetPtRatio.At(i))
                    lepton2jetP4 = self._p4_from_arrays(self.b_Muon_jetPt, self.b_Muon_jetEta, self.b_Muon_jetPhi, self.b_Muon_jetMass, i)
                    self.out_Lepton2_Jet_p4 = self._copy_p4_from_lv(lepton2jetP4)
                    self.out_Lepton2_mvaTOP = float(self.b_Muon_mvaTOP.At(i))
                    if self.isGoodMuon(i):
                        self.out_Lepton2_isTight = True
                    if self.isOnlyFakeMuon(i):
                        self.out_Lepton2_isFake = True
                if nNPMuons >= 2:
                    i = npMuonIdxs[1]
                    lepton3P4 = self._p4_from_arrays(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_phi, self.b_Muon_mass, i)
                    self.out_Lepton3_pdgId = -13 * int(self.b_Muon_charge.At(i))
                    self.out_Lepton3_p4 = self._copy_p4_from_lv(lepton3P4)
                    self.out_Lepton3_PFIso = float(self.b_Muon_pfRelIso04_all.At(i))
                    self.out_Lepton3_jetPtRatio = float(self.b_Muon_jetPtRatio.At(i))
                    lepton3jetP4 = self._p4_from_arrays(self.b_Muon_jetPt, self.b_Muon_jetEta, self.b_Muon_jetPhi, self.b_Muon_jetMass, i)
                    self.out_Lepton3_Jet_p4 = self._copy_p4_from_lv(lepton3jetP4)
                    self.out_Lepton3_mvaTOP = float(self.b_Muon_mvaTOP.At(i))
                    if self.isGoodMuon(i):
                        self.out_Lepton3_isTight = True
                    if self.isOnlyFakeMuon(i):
                        self.out_Lepton3_isFake = True

        self.out_GoodLeptonCode = 0
        if self.out_Lepton1_pdgId:
            self.out_GoodLeptonCode += 100
        if self.out_Lepton2_pdgId:
            self.out_GoodLeptonCode += 10
        if self.out_Lepton3_pdgId:
            self.out_GoodLeptonCode += 1
        if self.out_GoodLeptonCode != 111:
            return False

        if self.mode and self.mode.startswith("NPL"):
            if not (self.out_Lepton1_isFake or self.out_Lepton2_isFake or self.out_Lepton3_isFake):
                return False

        leptonTotP4 = lepton1P4 + lepton2P4 + lepton3P4
        self.out_TriLepton_mass = float(leptonTotP4.M())
        self.out_TriLepton_pt = float(leptonTotP4.Pt())

        if self.out_GoodLeptonCode >= 111 and self.mode in {"MuMuMu", "ElElEl", "TaTaTa", "NPLMuMuMu", "NPLElElEl", "NPLTaTaTa"}:
            c1x2 = self.out_Lepton1_pdgId * self.out_Lepton2_pdgId
            c2x3 = self.out_Lepton2_pdgId * self.out_Lepton3_pdgId
            c3x1 = self.out_Lepton3_pdgId * self.out_Lepton1_pdgId

            dm12 = abs((lepton1P4 + lepton2P4).M() - 91.2)
            dm23 = abs((lepton2P4 + lepton3P4).M() - 91.2)
            dm31 = abs((lepton3P4 + lepton1P4).M() - 91.2)

            select12 = 0
            select23 = 0
            select31 = 0

            if c1x2 < 0 and c2x3 > 0 and c3x1 > 0:
                select12 += 1
            elif c1x2 > 0 and c2x3 < 0 and c3x1 > 0:
                select23 += 1
            elif c1x2 > 0 and c2x3 > 0 and c3x1 < 0:
                select31 += 1
            elif c1x2 < 0 and c2x3 < 0 and c3x1 > 0:
                if dm12 < dm23:
                    select12 += 1
                else:
                    select23 += 1
            elif c1x2 > 0 and c2x3 < 0 and c3x1 < 0:
                if dm23 < dm31:
                    select23 += 1
                else:
                    select31 += 1
            elif c1x2 < 0 and c2x3 > 0 and c3x1 < 0:
                if dm12 < dm31:
                    select12 += 1
                else:
                    select31 += 1
            elif c1x2 < 0 and c2x3 < 0 and c3x1 < 0:
                if dm12 < dm23 and dm12 < dm31:
                    select12 += 1
                elif dm31 < dm12 and dm31 < dm23:
                    select31 += 1
                elif dm23 < dm12 and dm23 < dm31:
                    select23 += 1

            if select12 > 0:
                lepton1P4, lepton2P4, lepton3P4 = lepton3P4, lepton1P4, lepton2P4
                self._permute_lepton_aux_states((2, 0, 1))
            elif select31 > 0:
                lepton1P4, lepton2P4 = lepton2P4, lepton1P4
                self._permute_lepton_aux_states((1, 0, 2))

        if self.out_GoodLeptonCode >= 100:
            self.out_Lepton1_p4 = self._copy_p4_from_lv(lepton1P4)
        if self.out_GoodLeptonCode >= 110:
            self.out_Lepton2_p4 = self._copy_p4_from_lv(lepton2P4)
        if self.out_GoodLeptonCode >= 111:
            self.out_Lepton3_p4 = self._copy_p4_from_lv(lepton3P4)

        # Keep Lepton1/2/3 in their physics-role ordering:
        #   Lepton1 = W-associated lepton, Lepton2/3 = Z-candidate leptons.
        # Fill separate ptsortLepton1/2/3 branches in descending pT order.
        # These branches are used only for flavor/role-independent event-level
        # pT requirements and do not alter the original Lepton* branches.
        lepton_p4s = [lepton1P4, lepton2P4, lepton3P4]
        pt_order = sorted(
            range(3),
            key=lambda idx: float(lepton_p4s[idx].Pt()),
            reverse=True,
        )
        for rank, source_idx in enumerate(pt_order, start=1):
            source_p4 = lepton_p4s[source_idx]
            setattr(self, "out_ptsortLepton%d_pt" % rank, float(source_p4.Pt()))
            setattr(self, "out_ptsortLepton%d_eta" % rank, float(source_p4.Eta()))
            setattr(self, "out_ptsortLepton%d_phi" % rank, float(source_p4.Phi()))

        if self.out_Lepton2_pdgId * self.out_Lepton3_pdgId > 0:
            self.out_GoodLeptonCode *= -1

        if abs(self.out_GoodLeptonCode) % 100 == 11:
            zP4 = lepton2P4 + lepton3P4
            self.out_Z_p4 = self._copy_p4_from_lv(zP4)
            sumId = self.out_Lepton2_pdgId + self.out_Lepton3_pdgId
            if sumId <= (-11 - 11):
                self.out_Z_charge = +2
            elif sumId >= (+11 + 11):
                self.out_Z_charge = -2
            elif sumId <= -11:
                self.out_Z_charge = +1
            elif sumId >= +11:
                self.out_Z_charge = -1
            self.out_TriLepton_WleptonZdPhi = abs(lepton1P4.DeltaPhi(zP4))
            self.out_TriLepton_WleptonZdR = lepton1P4.DeltaR(zP4)

        if abs(self.out_GoodLeptonCode) >= 100:
            self.out_W_MT = self._compute_mt(lepton1P4, self.out_MET_pt, self.out_MET_phi)

        # Count loose UParTAK4 b-tagged jets with the AN2025/031 b-jet
        # acceptance: pT > 20 GeV and |eta| < 2.5.  These jets must pass the
        # same tight jet ID and lepton cleaning as the analysis jets, but the
        # b-tag pT threshold is intentionally lower than the GoodJet threshold.
        self.out_nBjet = 0
        for i in range(self.b_Jet_pt.GetSize()):
            if not self._is_btaggable_jet_for_tag_count(i):
                continue
            official_jet_id = self._passes_official_ak4puppi_tight_jet_id(i)
            if official_jet_id is not None:
                if not official_jet_id:
                    continue
            else:
                if not self._passes_tight_jet_id_from_composition(i):
                    continue
            jetP4ForBtag = self._p4_from_arrays(self.b_Jet_pt, self.b_Jet_eta, self.b_Jet_phi, self.b_Jet_mass, i)
            if lepton1P4.Pt() > 0 and lepton1P4.DeltaR(jetP4ForBtag) < 0.4:
                continue
            if lepton2P4.Pt() > 0 and lepton2P4.DeltaR(jetP4ForBtag) < 0.4:
                continue
            if lepton3P4.Pt() > 0 and lepton3P4.DeltaR(jetP4ForBtag) < 0.4:
                continue
            if self.b_Jet_btagDeepFlavB.At(i) > self.minBjetBDiscr_:
                self.out_nBjet += 1

        jetIdxs = []
        for i in range(self.b_Jet_pt.GetSize()):
            if not self.isGoodJet(i):
                continue
            jetP4 = self._p4_from_arrays(self.b_Jet_pt, self.b_Jet_eta, self.b_Jet_phi, self.b_Jet_mass, i)
            if lepton1P4.Pt() > 0 and lepton1P4.DeltaR(jetP4) < 0.4:
                continue
            if lepton2P4.Pt() > 0 and lepton2P4.DeltaR(jetP4) < 0.4:
                continue
            if lepton3P4.Pt() > 0 and lepton3P4.DeltaR(jetP4) < 0.4:
                continue
            jetIdxs.append(i)

            isTaggableJet = self._is_btaggable_jet_for_tag_count(i)
            isBTagged = isTaggableJet and self.b_Jet_btagDeepFlavB.At(i) > self.minBjetBDiscr_
            isCTagged = isTaggableJet and self.b_Jet_btagDeepFlavCvL.At(i) > self.minCjetCvLDiscr_ and self.b_Jet_btagDeepFlavCvB.At(i) > self.minCjetCvBDiscr_
            isQTagged = isTaggableJet and self.b_Jet_btagDeepFlavQG.At(i) > self.minQjetQGDiscr_

            if isBTagged:
                pass
            if isCTagged and not isBTagged:
                self.out_nCjet += 1
            if isQTagged and not isBTagged and not isCTagged:
                self.out_nUDSjet += 1

        self.out_nGoodJet = len(jetIdxs)
        jetIdxs.sort(key=lambda i: self.b_Jet_pt.At(i), reverse=True)
        for kk in jetIdxs:
            self.out_GoodJet_p4[0].append(float(self.b_Jet_pt.At(kk)))
            self.out_GoodJet_p4[1].append(float(self.b_Jet_eta.At(kk)))
            self.out_GoodJet_p4[2].append(float(self.b_Jet_phi.At(kk)))
            self.out_GoodJet_p4[3].append(float(self.b_Jet_mass.At(kk)))
            self.out_GoodJet_DeepFlavB.append(float(self.b_Jet_btagDeepFlavB.At(kk)))
            self.out_GoodJet_index.append(int(kk))
            self.out_GoodJet_DeepFlavCvL.append(float(self.b_Jet_btagDeepFlavCvL.At(kk)))
            self.out_GoodJet_DeepFlavCvB.append(float(self.b_Jet_btagDeepFlavCvB.At(kk)))
            self.out_GoodJet_DeepFlavQG.append(float(self.b_Jet_btagDeepFlavQG.At(kk)))
            if self.use_hadronFlavour or self._is_jetvar_mode():
                self.out_GoodJet_hadronFlavour.append(self._uchar_to_int(self.b_Jet_hadronFlavour.At(kk)))
            if self._is_jetvar_mode():
                self.out_GoodJet_Uncert.append(float(self.b_Jet_Uncert.At(kk)))

        tree = event._tree

        if not self._is_jetvar_mode():
            # Leptons
            for objName, p4, jetp4, pdgid, tau_decaymode, tau_prong, tau_genpartflav, pfiso, jetptratio, istight, isfake, jetb, mvatop in [
                ("Lepton1", self.out_Lepton1_p4, self.out_Lepton1_Jet_p4, self.out_Lepton1_pdgId, self.out_Lepton1_tauDecayMode, self.out_Lepton1_tauProng, self.out_Lepton1_tauGenPartFlav, self.out_Lepton1_PFIso, self.out_Lepton1_jetPtRatio, self.out_Lepton1_isTight, self.out_Lepton1_isFake, self.out_Lepton1_jetBTagDeepFlavor, self.out_Lepton1_mvaTOP),
                ("Lepton2", self.out_Lepton2_p4, self.out_Lepton2_Jet_p4, self.out_Lepton2_pdgId, self.out_Lepton2_tauDecayMode, self.out_Lepton2_tauProng, self.out_Lepton2_tauGenPartFlav, self.out_Lepton2_PFIso, self.out_Lepton2_jetPtRatio, self.out_Lepton2_isTight, self.out_Lepton2_isFake, self.out_Lepton2_jetBTagDeepFlavor, self.out_Lepton2_mvaTOP),
                ("Lepton3", self.out_Lepton3_p4, self.out_Lepton3_Jet_p4, self.out_Lepton3_pdgId, self.out_Lepton3_tauDecayMode, self.out_Lepton3_tauProng, self.out_Lepton3_tauGenPartFlav, self.out_Lepton3_PFIso, self.out_Lepton3_jetPtRatio, self.out_Lepton3_isTight, self.out_Lepton3_isFake, self.out_Lepton3_jetBTagDeepFlavor, self.out_Lepton3_mvaTOP),
            ]:
                out_map = {
                    "%s_pt" % objName: p4[0],
                    "%s_eta" % objName: p4[1],
                    "%s_phi" % objName: p4[2],
                    "%s_mass" % objName: p4[3],
                    "%s_PFIso" % objName: pfiso,
                    "%s_jetPtRatio" % objName: jetptratio,
                    "%s_Jet_pt" % objName: jetp4[0],
                    "%s_Jet_eta" % objName: jetp4[1],
                    "%s_Jet_phi" % objName: jetp4[2],
                    "%s_Jet_mass" % objName: jetp4[3],
                    "%s_jetBTagDeepFlavor" % objName: jetb,
                    "%s_pdgId" % objName: pdgid,
                    "%s_tauDecayMode" % objName: tau_decaymode,
                    "%s_tauProng" % objName: tau_prong,
                    "%s_tauGenPartFlav" % objName: tau_genpartflav,
                    "%s_isTight" % objName: istight,
                    "%s_isFake" % objName: isfake,
                    "%s_mvaTOP" % objName: mvatop,
                }
                for name, value in out_map.items():
                    self._set_tree_attr(tree, name, value)
                    self.out.fillBranch(name, value)

            ptsort_map = {
                "ptsortLepton1_pt": self.out_ptsortLepton1_pt,
                "ptsortLepton1_eta": self.out_ptsortLepton1_eta,
                "ptsortLepton1_phi": self.out_ptsortLepton1_phi,
                "ptsortLepton2_pt": self.out_ptsortLepton2_pt,
                "ptsortLepton2_eta": self.out_ptsortLepton2_eta,
                "ptsortLepton2_phi": self.out_ptsortLepton2_phi,
                "ptsortLepton3_pt": self.out_ptsortLepton3_pt,
                "ptsortLepton3_eta": self.out_ptsortLepton3_eta,
                "ptsortLepton3_phi": self.out_ptsortLepton3_phi,
            }
            for name, value in ptsort_map.items():
                self._set_tree_attr(tree, name, value)
                self.out.fillBranch(name, value)

            ptsort_goodjet_map = {
                "ptsortGoodJet1_pt": self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 1, 0),
                "ptsortGoodJet1_eta": self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 1, 1),
                "ptsortGoodJet1_phi": self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 1, 2),
                "ptsortGoodJet2_pt": self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 2, 0),
                "ptsortGoodJet2_eta": self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 2, 1),
                "ptsortGoodJet2_phi": self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 2, 2),
                "ptsortGoodJet3_pt": self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 3, 0),
                "ptsortGoodJet3_eta": self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 3, 1),
                "ptsortGoodJet3_phi": self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 3, 2),
            }
            for name, value in ptsort_goodjet_map.items():
                self._set_tree_attr(tree, name, value)
                self.out.fillBranch(name, value)

            other_p4_map = {
                "Z_pt": self.out_Z_p4[0],
                "Z_eta": self.out_Z_p4[1],
                "Z_phi": self.out_Z_p4[2],
                "Z_mass": self.out_Z_p4[3],
                "GoodJet_pt": self.out_GoodJet_p4[0],
                "GoodJet_eta": self.out_GoodJet_p4[1],
                "GoodJet_phi": self.out_GoodJet_p4[2],
                "GoodJet_mass": self.out_GoodJet_p4[3],
                "LeadingElectron_pt": self.out_LeadingElectron_p4[0],
                "LeadingElectron_eta": self.out_LeadingElectron_p4[1],
                "LeadingElectron_phi": self.out_LeadingElectron_p4[2],
                "LeadingElectron_mass": self.out_LeadingElectron_p4[3],
                "LeadingMuon_pt": self.out_LeadingMuon_p4[0],
                "LeadingMuon_eta": self.out_LeadingMuon_p4[1],
                "LeadingMuon_phi": self.out_LeadingMuon_p4[2],
                "LeadingMuon_mass": self.out_LeadingMuon_p4[3],
                "LeadingTau_pt": self.out_LeadingTau_p4[0],
                "LeadingTau_eta": self.out_LeadingTau_p4[1],
                "LeadingTau_phi": self.out_LeadingTau_p4[2],
                "LeadingTau_mass": self.out_LeadingTau_p4[3],
            }
            for name, value in other_p4_map.items():
                self._set_tree_attr(tree, name, value)
                self.out.fillBranch(name, value)

            scalar_map = {
                "MET_pt": self.out_MET_pt,
                "MET_phi": self.out_MET_phi,
                "LeadingLepton_pt": max(self.out_LeadingMuon_p4[0], self.out_LeadingElectron_p4[0], self.out_LeadingTau_p4[0]),
                "TriLepton_mass": self.out_TriLepton_mass,
                "TriLepton_pt": self.out_TriLepton_pt,
                "TriLepton_WleptonZdPhi": self.out_TriLepton_WleptonZdPhi,
                "TriLepton_WleptonZdR": self.out_TriLepton_WleptonZdR,
                "nVetoLepton": self.out_nVetoElectron + self.out_nVetoMuon,
                "nGoodElectron": self.out_nGoodElectron,
                "nGoodMuon": self.out_nGoodMuon,
                "nGoodTau": self.out_nGoodTau,
                "nNPTau": self.out_nNPTau,
                "nGoodLepton": self.out_nGoodElectron + self.out_nGoodMuon + self.out_nGoodTau,
                "nNPLepton": self.out_nNPElectron + self.out_nNPMuon + self.out_nNPTau,
                "GoodLeptonCode": self.out_GoodLeptonCode,
                "Z_charge": self.out_Z_charge,
                "W_MT": self.out_W_MT,
                "GoodJet_index": self.out_GoodJet_index,
                "GoodJet_DeepFlavB": self.out_GoodJet_DeepFlavB,
                **({"GoodJet_hadronFlavour": self.out_GoodJet_hadronFlavour} if self.use_hadronFlavour else {}),
                "nBjet": self.out_nBjet,
                "nCjet": self.out_nCjet,
                "nUDSjet": self.out_nUDSjet,
                "GoodJet_DeepFlavCvL": self.out_GoodJet_DeepFlavCvL,
                "GoodJet_DeepFlavCvB": self.out_GoodJet_DeepFlavCvB,
                "GoodJet_DeepFlavQG": self.out_GoodJet_DeepFlavQG,
            }
            for name, value in scalar_map.items():
                self._set_tree_attr(tree, name, value)
                self.out.fillBranch(name, value)

            self._set_tree_attr(tree, "nGoodJet", self.out_nGoodJet)
        else:
            suffix = self._jet_suffix()
            ptsort_goodjet_map = {
                "ptsortGoodJet1_pt_%s" % suffix: self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 1, 0),
                "ptsortGoodJet1_eta_%s" % suffix: self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 1, 1),
                "ptsortGoodJet1_phi_%s" % suffix: self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 1, 2),
                "ptsortGoodJet2_pt_%s" % suffix: self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 2, 0),
                "ptsortGoodJet2_eta_%s" % suffix: self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 2, 1),
                "ptsortGoodJet2_phi_%s" % suffix: self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 2, 2),
                "ptsortGoodJet3_pt_%s" % suffix: self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 3, 0),
                "ptsortGoodJet3_eta_%s" % suffix: self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 3, 1),
                "ptsortGoodJet3_phi_%s" % suffix: self._get_ptsort_goodjet(self.out_GoodJet_p4, self.out_nGoodJet, 3, 2),
            }
            for name, value in ptsort_goodjet_map.items():
                self._set_tree_attr(tree, name, value)
                self.out.fillBranch(name, value)

            other_map = {
                "GoodJet_pt_%s" % suffix: self._get_shifted_goodjet_pt(),
                "GoodJet_eta_%s" % suffix: self.out_GoodJet_p4[1],
                "GoodJet_phi_%s" % suffix: self.out_GoodJet_p4[2],
                "GoodJet_mass_%s" % suffix: self._get_shifted_goodjet_mass(),
                "W_MT_%s" % suffix: self.out_W_MT,
                "GoodJet_index_%s" % suffix: self.out_GoodJet_index,
                "GoodJet_DeepFlavB_%s" % suffix: self.out_GoodJet_DeepFlavB,
                "GoodJet_hadronFlavour_%s" % suffix: self.out_GoodJet_hadronFlavour,
                "nBjet_%s" % suffix: self.out_nBjet,
                "nCjet_%s" % suffix: self.out_nCjet,
                "nUDSjet_%s" % suffix: self.out_nUDSjet,
                "GoodJet_DeepFlavCvL_%s" % suffix: self.out_GoodJet_DeepFlavCvL,
                "GoodJet_DeepFlavCvB_%s" % suffix: self.out_GoodJet_DeepFlavCvB,
                "GoodJet_DeepFlavQG_%s" % suffix: self.out_GoodJet_DeepFlavQG,
            }
            for name, value in other_map.items():
                self._set_tree_attr(tree, name, value)
                self.out.fillBranch(name, value)

            corrected_met_pt, corrected_met_phi = self._get_corrected_met()
            self._set_tree_attr(tree, "MET_pt_%s" % suffix, corrected_met_pt)
            self.out.fillBranch("MET_pt_%s" % suffix, corrected_met_pt)
            self._set_tree_attr(tree, "MET_phi_%s" % suffix, corrected_met_phi)
            self.out.fillBranch("MET_phi_%s" % suffix, corrected_met_phi)
            self._set_tree_attr(tree, "nGoodJet_%s" % suffix, self.out_nGoodJet)
            self.out.fillBranch("nGoodJet_%s" % suffix, self.out_nGoodJet)

        return True

def _factory(mode, use_hadronFlavour=False, jetVar=None, jetSyst=None, year=None):
    year = str(year) if year is not None else "2024"
    electron_prompt_mva_cut = 0.64 if year == "2018" else 0.9
    muon_prompt_mva_cut = 0.64 if year == "2018" else 0.7

    # This module is used here for 2024 NanoAODv15.  For 2024, use the
    # UParTAK4 loose WP from the BTV JSON.  Keep the legacy 2018 value only for
    # backward-compatible factory names that may still be imported elsewhere.
    btag_loose_wp = 0.2783 if year == "2018" else FCNCTriLepton.UParTAK4_BTAG_WP_2024["L"]

    return lambda: FCNCTriLepton(
        mode=mode,
        use_hadronFlavour=use_hadronFlavour,
        jetVar=jetVar,
        jetSyst=jetSyst,
        muonPromptMVACut=muon_prompt_mva_cut,
        electronPromptMVACut=electron_prompt_mva_cut,
        minJetPt=30.0,
        maxJetAbsEta=4.7,
        applyPileupJetID=False,
        applyOfficialJetID=True,
        jetIDWorkingPoint="AK4PUPPI_Tight",
        btagWP="L",
        btagLooseWP=btag_loose_wp,
        maxBTagJetAbsEta=2.5,
    )


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

MODES = [
    "ElElEl", "MuElEl", "ElMuMu", "MuMuMu",
    "TaTaTa", "ElTaTa", "MuTaTa", "TaElEl", "TaMuMu",
    "NPLElElEl", "NPLMuElEl", "NPLElMuMu", "NPLMuMuMu",
    "NPLTaTaTa", "NPLElTaTa", "NPLMuTaTa", "NPLTaElEl", "NPLTaMuMu",
]

YEARS = ("2018", "2022", "2023", "2024")

for _year in YEARS:
    for _mode in MODES:
        globals()[f"fcnc_{_mode}_{_year}"] = _factory(_mode, use_hadronFlavour=False, year=_year)
        globals()[f"central_fcnc_{_mode}_{_year}"] = _factory(_mode, use_hadronFlavour=True, year=_year)
    for _jetVar in JET_VARS:
        for _mode in MODES:
            globals()[f"jetUncertUp_fcnc_{_mode}_{_year}_{_jetVar}"] = _factory(
                _mode, use_hadronFlavour=True, jetVar=_jetVar, jetSyst="Up", year=_year
            )
            globals()[f"jetUncertDown_fcnc_{_mode}_{_year}_{_jetVar}"] = _factory(
                _mode, use_hadronFlavour=True, jetVar=_jetVar, jetSyst="Down", year=_year
            )
