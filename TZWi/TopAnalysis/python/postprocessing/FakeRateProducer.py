import math
import os

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class _HistogramLookup2D(object):
    """Lookup helper for the 2024 signed-eta fake-rate maps.

    Confirmed map layout:
      X axis = lepton pT [GeV]
      Y axis = signed lepton eta
    """

    def __init__(self, histogram, name):
        if histogram is None:
            raise RuntimeError(
                "[FakeRateProducer] Histogram '%s' is not loaded." % name
            )
        if not histogram.InheritsFrom("TH2"):
            raise RuntimeError(
                "[FakeRateProducer] Object '%s' is not a TH2." % name
            )

        self.histogram = histogram
        self.name = name

    @staticmethod
    def _axis_bounds(axis):
        nbins = int(axis.GetNbins())
        return (
            float(axis.GetBinLowEdge(1)),
            float(axis.GetBinUpEdge(nbins)),
            nbins,
        )

    @staticmethod
    def _clamp_to_axis(value, axis):
        value = float(value)
        low, high, _ = _HistogramLookup2D._axis_bounds(axis)
        epsilon = 1.0e-6 * max(1.0, abs(high - low))

        if value < low:
            return low + epsilon
        if value >= high:
            return high - epsilon
        return value

    def find_bin(self, pt, eta):
        x_value = self._clamp_to_axis(pt, self.histogram.GetXaxis())
        y_value = self._clamp_to_axis(eta, self.histogram.GetYaxis())

        x_bin = int(self.histogram.GetXaxis().FindBin(x_value))
        y_bin = int(self.histogram.GetYaxis().FindBin(y_value))
        return x_bin, y_bin

    def get_by_bin(self, x_bin, y_bin):
        value = float(self.histogram.GetBinContent(x_bin, y_bin))
        if not math.isfinite(value):
            raise RuntimeError(
                "[FakeRateProducer] Non-finite value in '%s', bin=(%d,%d)."
                % (self.name, x_bin, y_bin)
            )
        return value


class _FakeRateMap2024(object):
    """Reader and validator for one 2024 signed-eta fake-rate ROOT file."""

    HISTOGRAM_NAMES = {
        "electron": (
            "h_fr_el",
            "h_fr_el_statUp",
            "h_fr_el_statDown",
        ),
        "muon": (
            "h_fr_mu",
            "h_fr_mu_statUp",
            "h_fr_mu_statDown",
        ),
    }

    def __init__(self, root_file, flavor, expected_btag_category=None):
        if flavor not in self.HISTOGRAM_NAMES:
            raise RuntimeError(
                "[FakeRateProducer] Unsupported flavor: %s" % flavor
            )

        self.root_file = root_file
        self.flavor = flavor
        self.expected_btag_category = expected_btag_category
        self.central = None
        self.stat_up = None
        self.stat_down = None
        self.central_lookup = None
        self.up_lookup = None
        self.down_lookup = None
        self.metadata = ""

    @staticmethod
    def _clone_histogram(root_file, name):
        obj = root_file.Get(name)
        if not obj:
            raise RuntimeError(
                "[FakeRateProducer] Histogram '%s' is missing." % name
            )
        if not obj.InheritsFrom("TH2"):
            raise RuntimeError(
                "[FakeRateProducer] Object '%s' is not a TH2." % name
            )

        clone = obj.Clone(name + "__clone")
        clone.SetDirectory(0)
        return clone

    @staticmethod
    def _same_axis(axis_a, axis_b, tolerance=1.0e-10):
        if int(axis_a.GetNbins()) != int(axis_b.GetNbins()):
            return False

        for index in range(1, int(axis_a.GetNbins()) + 2):
            edge_a = float(axis_a.GetBinLowEdge(index))
            edge_b = float(axis_b.GetBinLowEdge(index))
            if abs(edge_a - edge_b) > tolerance:
                return False
        return True

    @classmethod
    def _same_binning(cls, histogram_a, histogram_b):
        return (
            cls._same_axis(
                histogram_a.GetXaxis(), histogram_b.GetXaxis()
            )
            and cls._same_axis(
                histogram_a.GetYaxis(), histogram_b.GetYaxis()
            )
        )

    def _validate_axis_layout(self):
        x_axis = self.central.GetXaxis()
        y_axis = self.central.GetYaxis()

        x_low = float(x_axis.GetBinLowEdge(1))
        x_high = float(x_axis.GetBinUpEdge(x_axis.GetNbins()))
        y_low = float(y_axis.GetBinLowEdge(1))
        y_high = float(y_axis.GetBinUpEdge(y_axis.GetNbins()))

        valid = (
            x_low >= 0.0
            and x_high >= 100.0
            and y_low < 0.0
            and y_high > 0.0
            and abs(y_low) <= 5.0
            and abs(y_high) <= 5.0
        )
        if not valid:
            raise RuntimeError(
                "[FakeRateProducer] Unexpected axis layout in %s. "
                "Expected X=pT and Y=signed eta, obtained "
                "X=[%g,%g], Y=[%g,%g]."
                % (self.root_file, x_low, x_high, y_low, y_high)
            )

        print(
            "[FakeRateProducer]   axis layout : "
            "X=pT [%g,%g] (%d bins), "
            "Y=signed eta [%g,%g] (%d bins)"
            % (
                x_low,
                x_high,
                int(x_axis.GetNbins()),
                y_low,
                y_high,
                int(y_axis.GetNbins()),
            )
        )

    def _validate_bin_values(self):
        minimum = None
        maximum = None

        for x_bin in range(1, self.central.GetNbinsX() + 1):
            for y_bin in range(1, self.central.GetNbinsY() + 1):
                central = float(
                    self.central.GetBinContent(x_bin, y_bin)
                )
                up = float(
                    self.stat_up.GetBinContent(x_bin, y_bin)
                )
                down = float(
                    self.stat_down.GetBinContent(x_bin, y_bin)
                )

                values = (central, up, down)
                if not all(math.isfinite(value) for value in values):
                    raise RuntimeError(
                        "[FakeRateProducer] Non-finite FR value in %s, "
                        "bin=(%d,%d): central=%s, up=%s, down=%s"
                        % (
                            self.root_file,
                            x_bin,
                            y_bin,
                            central,
                            up,
                            down,
                        )
                    )

                if not (0.0 <= down <= central <= up < 1.0):
                    raise RuntimeError(
                        "[FakeRateProducer] Invalid FR ordering/range in %s, "
                        "bin=(%d,%d): expected "
                        "0 <= down <= central <= up < 1, obtained "
                        "central=%s, up=%s, down=%s"
                        % (
                            self.root_file,
                            x_bin,
                            y_bin,
                            central,
                            up,
                            down,
                        )
                    )

                minimum = (
                    central if minimum is None else min(minimum, central)
                )
                maximum = (
                    central if maximum is None else max(maximum, central)
                )

        print(
            "[FakeRateProducer]   central FR range: %.8g to %.8g"
            % (minimum, maximum)
        )

        if maximum >= 0.8:
            print(
                "[FakeRateProducer] WARNING: %s map contains FR >= 0.8. "
                "The corresponding f/(1-f) weight is large and can "
                "produce a statistically unstable nonprompt estimate."
                % self.flavor
            )

    def load(self):
        if not os.path.exists(self.root_file):
            raise RuntimeError(
                "[FakeRateProducer] Fake-rate ROOT file does not exist: %s"
                % self.root_file
            )

        root_file = ROOT.TFile.Open(self.root_file, "READ")
        if not root_file or root_file.IsZombie():
            raise RuntimeError(
                "[FakeRateProducer] Could not open fake-rate ROOT file: %s"
                % self.root_file
            )

        central_name, up_name, down_name = self.HISTOGRAM_NAMES[
            self.flavor
        ]

        try:
            self.central = self._clone_histogram(
                root_file, central_name
            )
            self.stat_up = self._clone_histogram(
                root_file, up_name
            )
            self.stat_down = self._clone_histogram(
                root_file, down_name
            )

            metadata = root_file.Get("metadata")
            if metadata:
                try:
                    self.metadata = str(metadata.GetTitle())
                except Exception:
                    self.metadata = ""
        finally:
            root_file.Close()

        if not self._same_binning(self.central, self.stat_up):
            raise RuntimeError(
                "[FakeRateProducer] central/statUp binning mismatch in %s"
                % self.root_file
            )
        if not self._same_binning(self.central, self.stat_down):
            raise RuntimeError(
                "[FakeRateProducer] central/statDown binning mismatch in %s"
                % self.root_file
            )

        if not self.metadata:
            raise RuntimeError(
                "[FakeRateProducer] Missing metadata in fake-rate map: %s"
                % self.root_file
            )
        if "etaMode=signed" not in self.metadata:
            raise RuntimeError(
                "[FakeRateProducer] The input map is not marked as signed-eta: "
                "%s" % self.metadata
            )
        if self.expected_btag_category is not None:
            expected_token = "btagCategory=%s" % self.expected_btag_category
            if expected_token not in self.metadata:
                raise RuntimeError(
                    "[FakeRateProducer] Wrong b-tag category map. Expected '%s' "
                    "in metadata from %s, obtained: %s"
                    % (expected_token, self.root_file, self.metadata)
                )

        self._validate_axis_layout()
        self._validate_bin_values()

        self.central_lookup = _HistogramLookup2D(
            self.central, central_name
        )
        self.up_lookup = _HistogramLookup2D(
            self.stat_up, up_name
        )
        self.down_lookup = _HistogramLookup2D(
            self.stat_down, down_name
        )

        print(
            "[FakeRateProducer] Loaded 2024 %s FR map from %s"
            % (self.flavor, self.root_file)
        )
        print(
            "[FakeRateProducer]   central/statUp/statDown: %s, %s, %s"
            % (central_name, up_name, down_name)
        )
        if self.metadata:
            print(
                "[FakeRateProducer]   metadata: %s" % self.metadata
            )

    def get(self, pt, eta):
        x_bin, y_bin = self.central_lookup.find_bin(pt, eta)

        central = self.central_lookup.get_by_bin(x_bin, y_bin)
        up = self.up_lookup.get_by_bin(x_bin, y_bin)
        down = self.down_lookup.get_by_bin(x_bin, y_bin)

        return {
            "flavor": self.flavor,
            "btag_category": self.expected_btag_category,
            "x_bin": x_bin,
            "y_bin": y_bin,
            "central": central,
            "up": up,
            "down": down,
            "err_up": abs(up - central),
            "err_down": abs(central - down),
        }


class FakeRateProducer(Module):
    """Apply either inclusive or b-tag-split 2024 signed-eta fake-rate maps.

    ``FAKERATE_MAP_MODE`` controls the map family:

      * ``inclusive``:
          one electron map and one muon map are used for every event.

      * ``split``:
          nBjet == 0 uses the 0B maps;
          nBjet >= 1 uses the GE1B maps.

    FCNCTriLepton must run before this producer when split mode is used,
    because split mode reads its cached ``b_out_nBjet`` value.

    Event weight:
      one fake   : + f/(1-f)
      two fakes  : - product_i[f_i/(1-f_i)]
      three fakes: + product_i[f_i/(1-f_i)]
    """

    MAP_MODE_INCLUSIVE = "inclusive"
    MAP_MODE_SPLIT = "split"

    CATEGORY_INCLUSIVE = "inclusive"
    CATEGORY_0B = "0b"
    CATEGORY_GE1B = "ge1b"

    def __init__(
        self,
        era=None,
        map_mode=None,
        fake_rate_file=None,
        electron_fake_rate_file=None,
        muon_fake_rate_file=None,
        electron_fake_rate_file_inclusive=None,
        muon_fake_rate_file_inclusive=None,
        electron_fake_rate_file_0b=None,
        electron_fake_rate_file_ge1b=None,
        muon_fake_rate_file_0b=None,
        muon_fake_rate_file_ge1b=None,
        use_abs_eta=False,
    ):
        self.era = (
            era
            if era is not None
            else os.environ.get("FAKERATE_ERA", "2024")
        )
        self.map_mode = self._normalize_map_mode(
            map_mode
            if map_mode is not None
            else os.environ.get("FAKERATE_MAP_MODE", self.MAP_MODE_SPLIT)
        )

        common_file = (
            fake_rate_file
            or os.environ.get("FAKERATE_2024_ROOT")
            or os.environ.get("FAKERATE_FILE")
        )
        common_el = (
            electron_fake_rate_file
            or os.environ.get("FAKERATE_ELECTRON_2024_ROOT")
            or os.environ.get("FAKERATE_ELECTRON_FILE")
            or common_file
        )
        common_mu = (
            muon_fake_rate_file
            or os.environ.get("FAKERATE_MUON_2024_ROOT")
            or os.environ.get("FAKERATE_MUON_FILE")
            or common_file
        )

        self.fake_rate_files = {
            self.CATEGORY_INCLUSIVE: {
                "electron": (
                    electron_fake_rate_file_inclusive
                    or os.environ.get("FAKERATE_ELECTRON_2024_INCLUSIVE_ROOT")
                    or os.environ.get("FAKERATE_ELECTRON_INCLUSIVE_ROOT")
                    or common_el
                ),
                "muon": (
                    muon_fake_rate_file_inclusive
                    or os.environ.get("FAKERATE_MUON_2024_INCLUSIVE_ROOT")
                    or os.environ.get("FAKERATE_MUON_INCLUSIVE_ROOT")
                    or common_mu
                ),
            },
            self.CATEGORY_0B: {
                "electron": (
                    electron_fake_rate_file_0b
                    or os.environ.get("FAKERATE_ELECTRON_2024_0B_ROOT")
                    or os.environ.get("FAKERATE_ELECTRON_0B_ROOT")
                    or common_el
                ),
                "muon": (
                    muon_fake_rate_file_0b
                    or os.environ.get("FAKERATE_MUON_2024_0B_ROOT")
                    or os.environ.get("FAKERATE_MUON_0B_ROOT")
                    or common_mu
                ),
            },
            self.CATEGORY_GE1B: {
                "electron": (
                    electron_fake_rate_file_ge1b
                    or os.environ.get("FAKERATE_ELECTRON_2024_GE1B_ROOT")
                    or os.environ.get("FAKERATE_ELECTRON_GE1B_ROOT")
                    or common_el
                ),
                "muon": (
                    muon_fake_rate_file_ge1b
                    or os.environ.get("FAKERATE_MUON_2024_GE1B_ROOT")
                    or os.environ.get("FAKERATE_MUON_GE1B_ROOT")
                    or common_mu
                ),
            },
        }

        requested_abs_eta = bool(
            use_abs_eta
            or os.environ.get(
                "FAKERATE_USE_ABS_ETA", ""
            ).lower() in ("1", "true", "yes")
        )
        if requested_abs_eta:
            raise RuntimeError(
                "[FakeRateProducer] The supplied 2024 maps use signed eta. "
                "FAKERATE_USE_ABS_ETA/use_abs_eta must not be enabled."
            )

        self._workers = {}

    @classmethod
    def _normalize_map_mode(cls, value):
        value = str(value or "").strip().lower()
        if value in (
            "split",
            "btag",
            "btag_split",
            "btag-split",
            "1",
            "true",
            "yes",
        ):
            return cls.MAP_MODE_SPLIT
        if value in (
            "inclusive",
            "inc",
            "0",
            "false",
            "no",
        ):
            return cls.MAP_MODE_INCLUSIVE
        raise RuntimeError(
            "[FakeRateProducer] Unsupported FAKERATE_MAP_MODE='%s'. "
            "Use split or inclusive." % value
        )

    @staticmethod
    def _topanalysis_data_dir(*parts):
        cmssw_base = os.environ.get("CMSSW_BASE", "")
        if not cmssw_base:
            return os.path.join(*parts)
        return os.path.join(
            cmssw_base,
            "src",
            "TZWi",
            "TopAnalysis",
            "data",
            *parts
        )

    @classmethod
    def _default_fake_rate_file(cls, flavor, category):
        if flavor not in ("electron", "muon"):
            raise RuntimeError(
                "[FakeRateProducer] Unsupported flavor: %s" % flavor
            )

        suffixes = {
            cls.CATEGORY_INCLUSIVE: "Inclusive",
            cls.CATEGORY_0B: "0B",
            cls.CATEGORY_GE1B: "GE1B",
        }
        if category not in suffixes:
            raise RuntimeError(
                "[FakeRateProducer] Unsupported map category: %s"
                % category
            )

        filename = (
            "fakerate_%s_promptMVA_2024_JetMET_SignedEta_%s.root"
            % (flavor, suffixes[category])
        )
        return cls._topanalysis_data_dir(
            "leptonFR",
            "2024",
            filename,
        )

    @classmethod
    def _resolve_fake_rate_file(cls, path, flavor, category):
        path = path or cls._default_fake_rate_file(flavor, category)

        if os.path.isabs(path):
            return path
        if os.path.exists(path):
            return os.path.abspath(path)

        candidates = (
            cls._topanalysis_data_dir(path),
            cls._topanalysis_data_dir("leptonFR", "2024", path),
            cls._topanalysis_data_dir("leptonFR", path),
        )
        for candidate in candidates:
            if os.path.exists(candidate):
                return candidate
        return path

    def _active_categories(self):
        if self.map_mode == self.MAP_MODE_INCLUSIVE:
            return (self.CATEGORY_INCLUSIVE,)
        return (self.CATEGORY_0B, self.CATEGORY_GE1B)

    def beginJob(self):
        if str(self.era) != "2024":
            raise RuntimeError(
                "[FakeRateProducer] This producer supports only era 2024; "
                "requested era=%s." % self.era
            )

        categories = self._active_categories()
        resolved = {}

        for category in categories:
            resolved[category] = {}
            for flavor in ("electron", "muon"):
                resolved[category][flavor] = self._resolve_fake_rate_file(
                    self.fake_rate_files[category][flavor],
                    flavor,
                    category,
                )

        if self.map_mode == self.MAP_MODE_SPLIT:
            allow_shared = os.environ.get(
                "FAKERATE_ALLOW_SHARED_BTAG_MAP",
                "0",
            ).lower() in ("1", "true", "yes")

            for flavor in ("electron", "muon"):
                file_0b = os.path.realpath(
                    resolved[self.CATEGORY_0B][flavor]
                )
                file_ge1b = os.path.realpath(
                    resolved[self.CATEGORY_GE1B][flavor]
                )
                if file_0b == file_ge1b and not allow_shared:
                    raise RuntimeError(
                        "[FakeRateProducer] The 0B and GE1B %s maps resolve "
                        "to the same file: %s." % (flavor, file_0b)
                    )

        self._workers = {}
        for category in categories:
            self._workers[category] = {}
            for flavor in ("electron", "muon"):
                worker = _FakeRateMap2024(
                    resolved[category][flavor],
                    flavor,
                    expected_btag_category=category,
                )
                worker.load()
                self._workers[category][flavor] = worker

        print(
            "[FakeRateProducer] FAKERATE_MAP_MODE=%s"
            % self.map_mode
        )
        if self.map_mode == self.MAP_MODE_INCLUSIVE:
            print(
                "[FakeRateProducer] Every event uses the Inclusive maps."
            )
        else:
            print("[FakeRateProducer] Event map selection:")
            print("[FakeRateProducer]   nBjet == 0 -> 0B maps")
            print("[FakeRateProducer]   nBjet >= 1 -> GE1B maps")

    def endJob(self):
        self._workers = {}

    def beginFile(
        self,
        inputFile,
        outputFile,
        inputTree,
        wrappedOutputTree,
    ):
        self.out = wrappedOutputTree

        for index in (1, 2, 3):
            self.out.branch("Lepton%d_FR" % index, "F")
            self.out.branch("Lepton%d_FRerrUp" % index, "F")
            self.out.branch("Lepton%d_FRerrDown" % index, "F")

        self.out.branch("FakeRateWeight", "F")
        self.out.branch("FakeRateWeightErrorUp", "F")
        self.out.branch("FakeRateWeightErrorDown", "F")

        # FakeRateMapMode:
        #   0 = inclusive maps
        #   1 = 0B/GE1B split maps
        self.out.branch("FakeRateMapMode", "I")

        # FakeRateBTagCategory:
        #  -1 = inclusive map
        #   0 = 0B map
        #   1 = GE1B map
        self.out.branch("FakeRateBTagCategory", "I")

    def endFile(
        self,
        inputFile,
        outputFile,
        inputTree,
        wrappedOutputTree,
    ):
        pass

    @staticmethod
    def _read_selected_lepton(tree, index):
        required = (
            "b_out_Lepton%d_pdgId" % index,
            "b_out_Lepton%d_pt" % index,
            "b_out_Lepton%d_eta" % index,
            "b_out_Lepton%d_isFake" % index,
        )
        missing = [
            name for name in required
            if not hasattr(tree, name)
        ]
        if missing:
            raise RuntimeError(
                "[FakeRateProducer] Missing cached selected-lepton "
                "value(s): %s. FCNCTriLepton must run before "
                "FakeRateProducer." % ", ".join(missing)
            )

        return (
            abs(int(getattr(tree, required[0]))),
            float(getattr(tree, required[1])),
            float(getattr(tree, required[2])),
            bool(getattr(tree, required[3])),
        )

    @classmethod
    def _read_btag_category(cls, tree):
        attr_name = "b_out_nBjet"
        if not hasattr(tree, attr_name):
            raise RuntimeError(
                "[FakeRateProducer] Missing cached nBjet value '%s'. "
                "Split mode requires FCNCTriLepton to run first."
                % attr_name
            )

        n_bjet = int(getattr(tree, attr_name))
        if n_bjet < 0:
            raise RuntimeError(
                "[FakeRateProducer] Invalid negative nBjet=%d."
                % n_bjet
            )

        category = (
            cls.CATEGORY_0B
            if n_bjet == 0
            else cls.CATEGORY_GE1B
        )
        return n_bjet, category

    def _select_event_map_category(self, tree):
        if self.map_mode == self.MAP_MODE_INCLUSIVE:
            n_bjet = (
                int(getattr(tree, "b_out_nBjet"))
                if hasattr(tree, "b_out_nBjet")
                else -1
            )
            return n_bjet, self.CATEGORY_INCLUSIVE

        return self._read_btag_category(tree)

    def _get_lepton_fr(self, pdgid, pt, eta, map_category):
        if map_category not in self._workers:
            raise RuntimeError(
                "[FakeRateProducer] No map set loaded for category %s."
                % map_category
            )

        if pdgid == 11:
            flavor = "electron"
        elif pdgid == 13:
            flavor = "muon"
        else:
            raise RuntimeError(
                "[FakeRateProducer] Unsupported selected-lepton "
                "pdgId=%s. This producer supports electrons and "
                "muons only." % pdgid
            )

        record = self._workers[map_category][flavor].get(pt, eta)
        record["btag_category"] = map_category
        return record

    @staticmethod
    def _validate_fr(fr):
        fr = float(fr)
        if not math.isfinite(fr) or not (0.0 <= fr < 1.0):
            raise RuntimeError(
                "[FakeRateProducer] Invalid fake rate: %s; "
                "expected 0 <= f < 1." % fr
            )
        return fr

    @classmethod
    def _odds(cls, fr):
        fr = cls._validate_fr(fr)
        return fr / (1.0 - fr)

    @classmethod
    def _event_weight(cls, fake_rates):
        number_of_fakes = len(fake_rates)
        if number_of_fakes == 0:
            return 0.0

        weight = 1.0
        for fr in fake_rates:
            weight *= cls._odds(fr)

        if number_of_fakes % 2 == 0:
            weight *= -1.0
        return float(weight)

    @classmethod
    def _exact_event_uncertainties(cls, fake_records):
        if not fake_records:
            return 0.0, 0.0, 0.0

        central_rates = [
            record["central"] for record in fake_records
        ]
        central_weight = cls._event_weight(central_rates)

        nuisance_groups = {}
        for index, record in enumerate(fake_records):
            key = (
                record["btag_category"],
                record["flavor"],
                record["x_bin"],
                record["y_bin"],
            )
            nuisance_groups.setdefault(key, []).append(index)

        error_up_squared = 0.0
        error_down_squared = 0.0

        for indices in nuisance_groups.values():
            up_rates = list(central_rates)
            down_rates = list(central_rates)

            for index in indices:
                up_rates[index] = fake_records[index]["up"]
                down_rates[index] = fake_records[index]["down"]

            up_weight = cls._event_weight(up_rates)
            down_weight = cls._event_weight(down_rates)

            error_up_squared += (
                up_weight - central_weight
            ) ** 2
            error_down_squared += (
                central_weight - down_weight
            ) ** 2

        return (
            central_weight,
            math.sqrt(error_up_squared),
            math.sqrt(error_down_squared),
        )

    def analyze(self, event):
        tree = event._tree
        _, map_category = self._select_event_map_category(tree)

        lepton_records = []
        fake_records = []

        for index in (1, 2, 3):
            pdgid, pt, eta, is_fake = self._read_selected_lepton(
                tree,
                index,
            )
            record = self._get_lepton_fr(
                pdgid,
                pt,
                eta,
                map_category,
            )
            lepton_records.append(record)

            if is_fake:
                fake_records.append(record)

        (
            fake_rate_weight,
            fake_rate_weight_error_up,
            fake_rate_weight_error_down,
        ) = self._exact_event_uncertainties(fake_records)

        for index, record in enumerate(
            lepton_records,
            start=1,
        ):
            self.out.fillBranch(
                "Lepton%d_FR" % index,
                record["central"],
            )
            self.out.fillBranch(
                "Lepton%d_FRerrUp" % index,
                record["err_up"],
            )
            self.out.fillBranch(
                "Lepton%d_FRerrDown" % index,
                record["err_down"],
            )

        self.out.fillBranch(
            "FakeRateWeight",
            fake_rate_weight,
        )
        self.out.fillBranch(
            "FakeRateWeightErrorUp",
            fake_rate_weight_error_up,
        )
        self.out.fillBranch(
            "FakeRateWeightErrorDown",
            fake_rate_weight_error_down,
        )
        self.out.fillBranch(
            "FakeRateMapMode",
            0
            if self.map_mode == self.MAP_MODE_INCLUSIVE
            else 1,
        )

        category_code = {
            self.CATEGORY_INCLUSIVE: -1,
            self.CATEGORY_0B: 0,
            self.CATEGORY_GE1B: 1,
        }[map_category]
        self.out.fillBranch(
            "FakeRateBTagCategory",
            category_code,
        )
        return True


FRProducer = lambda: FakeRateProducer()
FRProducer_2024 = lambda: FakeRateProducer(era="2024")
