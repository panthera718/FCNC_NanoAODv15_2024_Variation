import os
import math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class _HistogramLookup2D(object):
    """Safe 2D histogram lookup for fake-rate maps.

    The 2024 fake-rate maps produced from FR NTuples use
      X axis = lepton pT [GeV]
      Y axis = |eta|
    """

    def __init__(self, hist, name=""):
        self.hist = hist
        self.name = name

    @staticmethod
    def _clamp_to_axis(value, axis):
        value = float(value)
        nbins = axis.GetNbins()
        low = float(axis.GetBinLowEdge(1))
        high = float(axis.GetBinUpEdge(nbins))

        # Keep values inside the visible bins.  This avoids underflow/overflow
        # and makes values outside the measured range use the closest edge bin.
        eps = 1.0e-6 * max(1.0, abs(high - low))
        if value < low:
            return low + eps
        if value >= high:
            return high - eps
        return value

    def get(self, pt, abs_eta):
        if self.hist is None:
            print("[FakeRateProducer] ERROR: histogram '%s' is not loaded. Returning 0." % self.name)
            return 0.0

        x = self._clamp_to_axis(pt, self.hist.GetXaxis())
        y = self._clamp_to_axis(abs_eta, self.hist.GetYaxis())

        binx = self.hist.GetXaxis().FindBin(x)
        biny = self.hist.GetYaxis().FindBin(y)
        return float(self.hist.GetBinContent(binx, biny))


class _FakeRateMap2024(object):
    """Read central/up/down fake-rate maps from lepton_fake_rate_2024.root."""

    def __init__(self, root_file, channel):
        self.root_file = root_file
        self.channel = channel
        self.file_handle = None
        self.h_fr = None
        self.h_err_up = None
        self.h_err_down = None

    @staticmethod
    def _clone_object(tf, candidates, required=True):
        for name in candidates:
            obj = tf.Get(name)
            if obj:
                clone = obj.Clone(name.replace("/", "__") + "__clone")
                clone.SetDirectory(0)
                return clone
        if required:
            raise RuntimeError(
                "Could not find any of the following objects in fake-rate ROOT file: %s" %
                ", ".join(candidates)
            )
        return None

    def load(self):
        self.file_handle = ROOT.TFile.Open(self.root_file, "READ")
        if not self.file_handle or self.file_handle.IsZombie():
            raise RuntimeError("Could not open fake-rate ROOT file: %s" % self.root_file)

        if self.channel == "FREl":
            central_candidates = ["FREl/h_fr_FREl", "h_fr_FREl"]
            err_up_candidates = ["FREl/h_fr_errUp_FREl", "h_fr_errUp_FREl"]
            err_down_candidates = ["FREl/h_fr_errLow_FREl", "h_fr_errLow_FREl"]
        elif self.channel == "FRMu":
            central_candidates = ["FRMu/h_fr_FRMu", "h_fr_FRMu"]
            err_up_candidates = ["FRMu/h_fr_errUp_FRMu", "h_fr_errUp_FRMu"]
            err_down_candidates = ["FRMu/h_fr_errLow_FRMu", "h_fr_errLow_FRMu"]
        else:
            raise RuntimeError("Unsupported fake-rate channel: %s" % self.channel)

        self.h_fr = self._clone_object(self.file_handle, central_candidates, required=True)
        self.h_err_up = self._clone_object(self.file_handle, err_up_candidates, required=True)
        self.h_err_down = self._clone_object(self.file_handle, err_down_candidates, required=True)

        self.file_handle.Close()
        self.file_handle = None

        print("[FakeRateProducer] Loaded 2024 %s central FR map from %s" % (self.channel, self.root_file))
        print("[FakeRateProducer]   central hist : %s" % self.h_fr.GetName())
        print("[FakeRateProducer]   err-up hist  : %s" % self.h_err_up.GetName())
        print("[FakeRateProducer]   err-down hist: %s" % self.h_err_down.GetName())
        print("[FakeRateProducer]   lookup axes  : X = lepton pT, Y = |eta|")

    def get(self, pt, abs_eta):
        fr = _HistogramLookup2D(self.h_fr, self.channel + " central").get(pt, abs_eta)
        err_up = _HistogramLookup2D(self.h_err_up, self.channel + " errUp").get(pt, abs_eta)
        err_down = _HistogramLookup2D(self.h_err_down, self.channel + " errDown").get(pt, abs_eta)
        return float(fr), float(err_up), float(err_down)


class FakeRateProducer(Module, object):
    """Apply 2024 electron/muon fake rates to NPL events.

    The expected input fake-rate ROOT file is the output of the FR measurement
    workflow, with objects such as:
      FREl/h_fr_FREl, FREl/h_fr_errUp_FREl, FREl/h_fr_errLow_FREl
      FRMu/h_fr_FRMu, FRMu/h_fr_errUp_FRMu, FRMu/h_fr_errLow_FRMu

    The lookup variables are:
      electron: Lepton{i}_pt, abs(Lepton{i}_eta)
      muon    : Lepton{i}_pt, abs(Lepton{i}_eta)

    This intentionally does not use pTcone because the newly measured 2024 maps
    were filled with FR_lepton_pt and FR_lepton_absEta.
    """

    def __init__(self, era=None, fake_rate_file=None):
        self.era = era if era is not None else os.environ.get("FAKERATE_ERA", "2024")
        self.fake_rate_file = fake_rate_file or os.environ.get("FAKERATE_2024_ROOT") or os.environ.get("FAKERATE_FILE")
        self._worker_el = None
        self._worker_mu = None
        self.inputFileName = ""

    @staticmethod
    def _default_fake_rate_file():
        cmssw_base = os.environ.get("CMSSW_BASE", "")
        if not cmssw_base:
            return "lepton_fake_rate_2024.root"
        return os.path.join(
            cmssw_base,
            "src/TZWi/TopAnalysis/data/leptonFR/lepton_fake_rate_2024.root",
        )

    def _resolve_fake_rate_file(self):
        path = self.fake_rate_file or self._default_fake_rate_file()
        if os.path.isabs(path):
            return path

        # First try the path as given relative to the current working directory.
        if os.path.exists(path):
            return os.path.abspath(path)

        # Then try the standard TopAnalysis data/leptonFR directory.
        cmssw_base = os.environ.get("CMSSW_BASE", "")
        if cmssw_base:
            candidate = os.path.join(cmssw_base, "src/TZWi/TopAnalysis/data/leptonFR", path)
            if os.path.exists(candidate):
                return candidate

        # Return the original path so the error message shows what was requested.
        return path

    def beginJob(self):
        if str(self.era) != "2024":
            print(
                "[FakeRateProducer] Warning: this rewritten producer is intended for 2024 "
                "electron/muon fake-rate maps. Requested era='%s'." % str(self.era)
            )

        fr_file = self._resolve_fake_rate_file()
        self._worker_el = _FakeRateMap2024(fr_file, "FREl")
        self._worker_mu = _FakeRateMap2024(fr_file, "FRMu")
        self._worker_el.load()
        self._worker_mu.load()

    def endJob(self):
        self._worker_el = None
        self._worker_mu = None

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.inputFileName = inputFile.GetName() if inputFile and hasattr(inputFile, "GetName") else ""

        # Keep central branch names consistent with the previous producer.
        self.out.branch("Lepton1_FR", "F")
        self.out.branch("Lepton2_FR", "F")
        self.out.branch("Lepton3_FR", "F")
        self.out.branch("FakeRateWeight", "F")

        # Replace the old symmetric error branches by asymmetric uncertainties.
        self.out.branch("Lepton1_FRerrUp", "F")
        self.out.branch("Lepton1_FRerrDown", "F")
        self.out.branch("Lepton2_FRerrUp", "F")
        self.out.branch("Lepton2_FRerrDown", "F")
        self.out.branch("Lepton3_FRerrUp", "F")
        self.out.branch("Lepton3_FRerrDown", "F")
        self.out.branch("FakeRateWeightErrorUp", "F")
        self.out.branch("FakeRateWeightErrorDown", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    @staticmethod
    def _get_lepton_pt_eta_from_tree(tree, idx):
        pt = float(getattr(tree, "b_out_Lepton%d_pt" % idx))
        abs_eta = abs(float(getattr(tree, "b_out_Lepton%d_eta" % idx)))
        return pt, abs_eta

    def _get_lepton_fr(self, tree, idx):
        pdgid = abs(int(getattr(tree, "b_out_Lepton%d_pdgId" % idx)))
        pt, abs_eta = self._get_lepton_pt_eta_from_tree(tree, idx)

        if pdgid == 11:
            return self._worker_el.get(pt, abs_eta)

        if pdgid == 13:
            return self._worker_mu.get(pt, abs_eta)

        # This producer was rewritten for electron/muon fake rates.  Tau fake
        # rates should be handled by a separate tau-specific producer if needed.
        return 0.0, 0.0, 0.0

    @staticmethod
    def _safe_fr(fr):
        fr = float(fr)
        if fr >= 1.0:
            return 1.0 - 1.0e-6
        if fr <= -1.0:
            return -1.0 + 1.0e-6
        return fr

    @classmethod
    def _fr_to_factor(cls, fr):
        fr = cls._safe_fr(fr)
        return -fr / (1.0 - fr)

    @classmethod
    def _frerr_to_factorerr(cls, fr, fr_err):
        fr = cls._safe_fr(fr)
        fr_err = abs(float(fr_err))
        return fr_err / ((1.0 - fr) ** 2)

    @staticmethod
    def _weight_error_from_terms(fake_terms, use_up=True):
        if not fake_terms:
            return 0.0

        err2 = 0.0
        for i, (_, factor_err_up, factor_err_down) in enumerate(fake_terms):
            factor_err = factor_err_up if use_up else factor_err_down
            other_prod = -1.0
            for j, (other_factor, _, _) in enumerate(fake_terms):
                if i == j:
                    continue
                other_prod *= other_factor
            err2 += (other_prod * factor_err) ** 2
        return math.sqrt(err2)

    def analyze(self, event):
        tree = event._tree

        fr1, fr1_up, fr1_down = self._get_lepton_fr(tree, 1)
        fr2, fr2_up, fr2_down = self._get_lepton_fr(tree, 2)
        fr3, fr3_up, fr3_down = self._get_lepton_fr(tree, 3)

        fake_terms = []
        for is_fake, fr, err_up, err_down in [
            (tree.b_out_Lepton1_isFake, fr1, fr1_up, fr1_down),
            (tree.b_out_Lepton2_isFake, fr2, fr2_up, fr2_down),
            (tree.b_out_Lepton3_isFake, fr3, fr3_up, fr3_down),
        ]:
            if is_fake:
                fake_terms.append((
                    self._fr_to_factor(fr),
                    self._frerr_to_factorerr(fr, err_up),
                    self._frerr_to_factorerr(fr, err_down),
                ))

        if not fake_terms:
            fr_weight = 0.0
            fr_weight_err_up = 0.0
            fr_weight_err_down = 0.0
        else:
            fr_weight = -1.0
            for factor, _, _ in fake_terms:
                fr_weight *= factor
            fr_weight_err_up = self._weight_error_from_terms(fake_terms, use_up=True)
            fr_weight_err_down = self._weight_error_from_terms(fake_terms, use_up=False)

        self.out.fillBranch("Lepton1_FR", fr1)
        self.out.fillBranch("Lepton1_FRerrUp", fr1_up)
        self.out.fillBranch("Lepton1_FRerrDown", fr1_down)
        self.out.fillBranch("Lepton2_FR", fr2)
        self.out.fillBranch("Lepton2_FRerrUp", fr2_up)
        self.out.fillBranch("Lepton2_FRerrDown", fr2_down)
        self.out.fillBranch("Lepton3_FR", fr3)
        self.out.fillBranch("Lepton3_FRerrUp", fr3_up)
        self.out.fillBranch("Lepton3_FRerrDown", fr3_down)
        self.out.fillBranch("FakeRateWeight", fr_weight)
        self.out.fillBranch("FakeRateWeightErrorUp", fr_weight_err_up)
        self.out.fillBranch("FakeRateWeightErrorDown", fr_weight_err_down)

        return True


FRProducer = lambda: FakeRateProducer()
FRProducer_2024 = lambda: FakeRateProducer(era="2024")
