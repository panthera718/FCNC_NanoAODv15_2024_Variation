import math
import os
from math import pi

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class _HistogramLookup2D(object):
    """
    Lookup helper for the supplied signed-eta promptMVA SF maps.

    Confirmed map layout:
      x-axis = lepton pT
      y-axis = signed eta
    """

    def __init__(self, histogram, label=""):
        self.histogram = histogram
        self.label = label

    @staticmethod
    def _clamp_to_axis(value, axis):
        value = float(value)
        nbins = int(axis.GetNbins())
        low = float(axis.GetBinLowEdge(1))
        high = float(axis.GetBinUpEdge(nbins))
        epsilon = 1.0e-6 * max(1.0, abs(high - low))

        if value < low:
            return low + epsilon
        if value >= high:
            return high - epsilon
        return value

    def _find_bin(self, pt, eta):
        if self.histogram is None:
            raise RuntimeError(
                "[lepSFProducer] Input SF histogram is unavailable: %s"
                % self.label
            )

        x_value = self._clamp_to_axis(
            pt, self.histogram.GetXaxis()
        )
        y_value = self._clamp_to_axis(
            eta, self.histogram.GetYaxis()
        )

        x_bin = int(self.histogram.GetXaxis().FindBin(x_value))
        y_bin = int(self.histogram.GetYaxis().FindBin(y_value))
        return x_bin, y_bin

    def get(self, pt, eta):
        x_bin, y_bin = self._find_bin(pt, eta)
        value = float(self.histogram.GetBinContent(x_bin, y_bin))

        if not math.isfinite(value):
            raise RuntimeError(
                "[lepSFProducer] Non-finite SF in %s at pT=%s, eta=%s, "
                "bin=(%d,%d)"
                % (self.label, pt, eta, x_bin, y_bin)
            )
#        if value <= 0.0:
        if value < 0.0:
            raise RuntimeError(
                "[lepSFProducer] Non-positive SF in %s at pT=%s, eta=%s, "
                "bin=(%d,%d): %s"
                % (self.label, pt, eta, x_bin, y_bin, value)
            )

        return value

    def get_error(self, pt, eta):
        x_bin, y_bin = self._find_bin(pt, eta)
        value = float(self.histogram.GetBinError(x_bin, y_bin))
        if not math.isfinite(value):
            raise RuntimeError(
                "[lepSFProducer] Non-finite SF error in %s at pT=%s, eta=%s"
                % (self.label, pt, eta)
            )
        return abs(value)


class _PromptMVALeptonSFMap(object):
    """
    Reader for the user-measured promptMVA lepton SF ROOT maps.

    Expected objects:
      h2_sf
      h2_sf_statUp
      h2_sf_statDown

    The supplied files have:
      x = pT
      y = signed eta
    """

    def __init__(self, root_file_name, hist_name="h2_sf", label=""):
        self.root_file_name = root_file_name
        self.hist_name = hist_name
        self.label = label if label else hist_name

        self.hist = None
        self.hist_stat_up = None
        self.hist_stat_down = None
        self.metadata = ""

        self._load()

    @staticmethod
    def _clone(root_file, name, required=True):
        obj = root_file.Get(name)
        if not obj:
            if required:
                raise RuntimeError(
                    "[lepSFProducer] Histogram %s is missing." % name
                )
            return None

        clone = obj.Clone(name.replace("/", "__") + "__clone")
        clone.SetDirectory(0)
        return clone

    @staticmethod
    def _axis_bounds(axis):
        nbins = int(axis.GetNbins())
        return (
            float(axis.GetBinLowEdge(1)),
            float(axis.GetBinUpEdge(nbins)),
            nbins,
        )

    @staticmethod
    def _histogram_range(histogram):
        minimum = None
        maximum = None

        for x_bin in range(1, histogram.GetNbinsX() + 1):
            for y_bin in range(1, histogram.GetNbinsY() + 1):
                value = float(histogram.GetBinContent(x_bin, y_bin))
                if not math.isfinite(value):
                    raise RuntimeError(
                        "[lepSFProducer] Non-finite value in %s, bin=(%d,%d)"
                        % (histogram.GetName(), x_bin, y_bin)
                    )
#                if value <= 0.0:
                if value < 0.0:
                    raise RuntimeError(
                        "[lepSFProducer] Non-positive value in %s, bin=(%d,%d): %s"
                        % (histogram.GetName(), x_bin, y_bin, value)
                    )
                minimum = value if minimum is None else min(minimum, value)
                maximum = value if maximum is None else max(maximum, value)

        return float(minimum), float(maximum)

    def _validate_axis_layout(self):
        x_low, x_high, x_nbins = self._axis_bounds(
            self.hist.GetXaxis()
        )
        y_low, y_high, y_nbins = self._axis_bounds(
            self.hist.GetYaxis()
        )

        expected_layout = (
            x_low >= 0.0
            and x_high >= 100.0
            and y_low < 0.0
            and y_high > 0.0
            and abs(y_low) <= 5.0
            and abs(y_high) <= 5.0
        )

        if not expected_layout:
            raise RuntimeError(
                "[lepSFProducer] Unexpected axis layout in %s/%s. "
                "Expected x=pT and y=signed eta, but obtained "
                "x=[%s,%s] (%d bins), y=[%s,%s] (%d bins)."
                % (
                    self.root_file_name,
                    self.hist_name,
                    x_low,
                    x_high,
                    x_nbins,
                    y_low,
                    y_high,
                    y_nbins,
                )
            )

        print(
            "[lepSFProducer] %s axis layout: x=pT [%g,%g] (%d bins), "
            "y=signed eta [%g,%g] (%d bins)"
            % (
                self.label,
                x_low,
                x_high,
                x_nbins,
                y_low,
                y_high,
                y_nbins,
            )
        )

    def _load(self):
        root_file = ROOT.TFile.Open(self.root_file_name, "READ")
        if not root_file or root_file.IsZombie():
            raise RuntimeError(
                "[lepSFProducer] Cannot open lepton SF ROOT file: %s"
                % self.root_file_name
            )

        try:
            self.hist = self._clone(
                root_file, self.hist_name, required=True
            )
            self.hist_stat_up = self._clone(
                root_file,
                self.hist_name + "_statUp",
                required=False,
            )
            self.hist_stat_down = self._clone(
                root_file,
                self.hist_name + "_statDown",
                required=False,
            )

            metadata = root_file.Get("metadata")
            if metadata:
                try:
                    self.metadata = str(metadata.GetTitle())
                except Exception:
                    self.metadata = ""
        finally:
            root_file.Close()

        for name, histogram in (
            (self.hist_name, self.hist),
            (self.hist_name + "_statUp", self.hist_stat_up),
            (self.hist_name + "_statDown", self.hist_stat_down),
        ):
            if histogram is not None and not histogram.InheritsFrom("TH2"):
                raise RuntimeError(
                    "[lepSFProducer] %s in %s is not a TH2."
                    % (name, self.root_file_name)
                )

        self._validate_axis_layout()

        central_min, central_max = self._histogram_range(self.hist)
        print(
            "[lepSFProducer] Loaded %s from %s"
            % (self.hist_name, self.root_file_name)
        )
        print(
            "[lepSFProducer]   central SF range: %.8g to %.8g"
            % (central_min, central_max)
        )
        print(
            "[lepSFProducer]   statUp: %s; statDown: %s"
            % (
                "loaded" if self.hist_stat_up else "using central bin errors",
                "loaded" if self.hist_stat_down else "using central bin errors",
            )
        )
        if self.metadata:
            print("[lepSFProducer]   metadata: %s" % self.metadata)

        # This is a diagnostic only.  The supplied maps are still applied as-is.
        if central_min < 0.5 or central_max > 1.5:
            print(
                "[lepSFProducer] WARNING: %s map contains values far from 1 "
                "(range %.6g to %.6g). This producer will apply the map "
                "exactly as supplied; verify the SF measurement if such "
                "values are not expected."
                % (self.label, central_min, central_max)
            )

    def get_triplet(self, pt, eta):
        central_lookup = _HistogramLookup2D(
            self.hist, self.label + " central"
        )
        central = central_lookup.get(pt, eta)

        if self.hist_stat_up is not None:
            up = _HistogramLookup2D(
                self.hist_stat_up, self.label + " statUp"
            ).get(pt, eta)
        else:
            up = central + central_lookup.get_error(pt, eta)

        if self.hist_stat_down is not None:
            down = _HistogramLookup2D(
                self.hist_stat_down, self.label + " statDown"
            ).get(pt, eta)
        else:
            down = central - central_lookup.get_error(pt, eta)

#        if down <= 0.0:
        if down < 0.0:
            raise RuntimeError(
                "[lepSFProducer] Non-positive down variation for %s at "
                "pT=%s, eta=%s: central=%s, down=%s"
                % (self.label, pt, eta, central, down)
            )

        error_up = abs(up - central)
        error_down = abs(central - down)

        return (
            float(central),
            float(up),
            float(down),
            float(error_up),
            float(error_down),
        )


class lepSFProducer(Module):
    """
    Apply the user-measured 2024 promptMVA/selection SF maps.

    Files:
      data/leptonSF/2024/leptonSF_electron_2024_signedEta.root
      data/leptonSF/2024/leptonSF_muon_2024_signedEta.root

    This module does not apply official electron reconstruction SFs.  Those are
    handled separately by lepRecoSFProducer.py.
    """

    def __init__(
        self,
        electron_file=None,
        muon_file=None,
        hist_name="h2_sf",
        maxMatchDR=0.02,
        maxMatchRelPt=0.02,
        maxFallbackWarnings=5,
    ):
        cmssw_base = os.environ.get("CMSSW_BASE", "")
        base = os.path.join(
            cmssw_base,
            "src", "TZWi", "TopAnalysis", "data", "leptonSF", "2024",
        )

        self.electron_file = (
            electron_file
            if electron_file is not None
            else os.path.join(
                base, "leptonSF_electron_2024_signedEta.root"
            )
        )
        self.muon_file = (
            muon_file
            if muon_file is not None
            else os.path.join(
                base, "leptonSF_muon_2024_signedEta.root"
            )
        )
        self.hist_name = hist_name

        self.max_match_dr = float(maxMatchDR)
        self.max_match_relpt = float(maxMatchRelPt)
        self.max_fallback_warnings = int(maxFallbackWarnings)
        self._fallback_warnings = 0

        self._worker_el = None
        self._worker_mu = None

    def beginJob(self):
        self._worker_el = _PromptMVALeptonSFMap(
            self.electron_file,
            self.hist_name,
            "electron",
        )
        self._worker_mu = _PromptMVALeptonSFMap(
            self.muon_file,
            self.hist_name,
            "muon",
        )

    def endJob(self):
        self._worker_el = None
        self._worker_mu = None

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        for prefix in ("Lepton", "Electron", "Muon"):
            self.out.branch(prefix + "_SF", "F")
            self.out.branch(prefix + "_SFUp", "F")
            self.out.branch(prefix + "_SFDown", "F")
            self.out.branch(prefix + "_SFerr", "F")
            self.out.branch(prefix + "_SFerrUp", "F")
            self.out.branch(prefix + "_SFerrDown", "F")
            self.out.branch(prefix + "_SFerrStat", "F")
            self.out.branch(prefix + "_SFerrSyst", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    @staticmethod
    def _delta_phi(phi1, phi2):
        dphi = float(phi1) - float(phi2)
        while dphi > pi:
            dphi -= 2.0 * pi
        while dphi <= -pi:
            dphi += 2.0 * pi
        return dphi

    @classmethod
    def _delta_r2(cls, eta1, phi1, eta2, phi2):
        deta = float(eta1) - float(eta2)
        dphi = cls._delta_phi(phi1, phi2)
        return deta * deta + dphi * dphi

    @staticmethod
    def _electron_sc_eta(electron):
        try:
            return float(electron.superclusterEta)
        except Exception:
            pass

        try:
            return float(electron.eta) + float(electron.deltaEtaSC)
        except Exception:
            return float(electron.eta)

    def _matched_electron_sc_eta(self, event, pdgid, pt, eta, phi):
        fallback_eta = float(eta)

        try:
            electrons = Collection(event, "Electron")
        except Exception:
            return fallback_eta

        target_pdgid = int(pdgid)
        best = None
        best_dr2 = None
        best_relpt = None
        best_score = None

        for electron in electrons:
            try:
                electron_pdgid = -11 * int(electron.charge)
                if electron_pdgid != target_pdgid:
                    continue

                dr2 = self._delta_r2(
                    eta, phi, electron.eta, electron.phi
                )
                relpt = abs(float(pt) - float(electron.pt)) / max(
                    1.0, abs(float(pt))
                )
                score = dr2 + relpt * relpt

                if best_score is None or score < best_score:
                    best = electron
                    best_dr2 = dr2
                    best_relpt = relpt
                    best_score = score
            except Exception:
                continue

        matched = (
            best is not None
            and best_dr2 is not None
            and best_relpt is not None
            and best_dr2 <= self.max_match_dr * self.max_match_dr
            and best_relpt <= self.max_match_relpt
        )

        if matched:
            return self._electron_sc_eta(best)

        if self._fallback_warnings < self.max_fallback_warnings:
            print(
                "[lepSFProducer] WARNING: no Electron collection match for "
                "selected electron pdgId=%d, pt=%.6g, eta=%.6g, phi=%.6g; "
                "using Lepton_eta for the promptMVA SF lookup."
                % (int(pdgid), float(pt), float(eta), float(phi))
            )
            self._fallback_warnings += 1

        return fallback_eta

    @staticmethod
    def _read_selected_lepton(event, index):
        tree = event._tree
        return (
            int(getattr(tree, "b_out_Lepton%d_pdgId" % index)),
            float(getattr(tree, "b_out_Lepton%d_pt" % index)),
            float(getattr(tree, "b_out_Lepton%d_eta" % index)),
            float(getattr(tree, "b_out_Lepton%d_phi" % index)),
        )

    @staticmethod
    def _combine_shifted_products(components):
        central = 1.0
        up = 1.0
        down = 1.0

        for sf, sf_up, sf_down in components:
            central *= float(sf)
            up *= float(sf_up)
            down *= float(sf_down)

        error_up = abs(up - central)
        error_down = abs(central - down)
        error = max(error_up, error_down)

        return (
            float(central),
            float(up),
            float(down),
            float(error),
            float(error_up),
            float(error_down),
        )

    def _process_one_lepton(self, event, index):
        pdgid, pt, eta, phi = self._read_selected_lepton(event, index)
        abs_pdgid = abs(pdgid)

        if abs_pdgid == 11:
            lookup_eta = self._matched_electron_sc_eta(
                event, pdgid, pt, eta, phi
            )
            sf, sf_up, sf_down, _, _ = self._worker_el.get_triplet(
                pt, lookup_eta
            )
            return "electron", sf, sf_up, sf_down

        if abs_pdgid == 13:
            sf, sf_up, sf_down, _, _ = self._worker_mu.get_triplet(
                pt, eta
            )
            return "muon", sf, sf_up, sf_down

        return None, 1.0, 1.0, 1.0

    def _fill_prefix(self, prefix, values):
        sf, sf_up, sf_down, error, error_up, error_down = values

        self.out.fillBranch(prefix + "_SF", sf)
        self.out.fillBranch(prefix + "_SFUp", sf_up)
        self.out.fillBranch(prefix + "_SFDown", sf_down)
        self.out.fillBranch(prefix + "_SFerr", error)
        self.out.fillBranch(prefix + "_SFerrUp", error_up)
        self.out.fillBranch(prefix + "_SFerrDown", error_down)
        self.out.fillBranch(prefix + "_SFerrStat", error)
        self.out.fillBranch(prefix + "_SFerrSyst", 0.0)

    def analyze(self, event):
        all_components = []
        electron_components = []
        muon_components = []

        for index in (1, 2, 3):
            flavor, sf, sf_up, sf_down = self._process_one_lepton(
                event, index
            )
            component = (sf, sf_up, sf_down)

            all_components.append(component)
            if flavor == "electron":
                electron_components.append(component)
            elif flavor == "muon":
                muon_components.append(component)

        self._fill_prefix(
            "Lepton",
            self._combine_shifted_products(all_components),
        )
        self._fill_prefix(
            "Electron",
            self._combine_shifted_products(electron_components),
        )
        self._fill_prefix(
            "Muon",
            self._combine_shifted_products(muon_components),
        )

        return True


lepSF = lambda: lepSFProducer()
lepSF_2024 = lambda: lepSFProducer()
