import bisect
import gzip
import json
import math
import os
from math import pi

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class _ElectronRecoJSONEvaluator(object):
    """
    Minimal evaluator for the official 2024 EGM Electron-ID-SF correction.

    The supplied electron.json.gz uses correctionlib schema v2 with

      year -> ValType -> WorkingPoint -> multibinning(eta, pt)

    and contains string edges "-inf" and "inf".  The correctionlib version
    shipped with CMSSW_13_3_3 can reject these valid edges with

      "binning edges are not monotone increasing"

    so this evaluator reads this specific official schema directly with the
    Python standard library.

    Multibinning flattening follows correctionlib schema ordering:
    the last input axis varies fastest.  For inputs ["eta", "pt"],

      flat_index = eta_bin * n_pt_bins + pt_bin
    """

    def __init__(self, json_path, correction_name="Electron-ID-SF"):
        self.json_path = json_path
        self.correction_name = correction_name
        self.correction = None
        self.input_names = None
        self._load()

    @staticmethod
    def _open_json(path):
        if path.endswith(".gz"):
            with gzip.open(path, "rt") as handle:
                return json.load(handle)
        with open(path, "rt") as handle:
            return json.load(handle)

    @staticmethod
    def _edge_to_float(edge):
        if edge == "-inf":
            return float("-inf")
        if edge == "inf":
            return float("inf")
        return float(edge)

    @classmethod
    def _validate_edges(cls, edges, label):
        converted = [cls._edge_to_float(edge) for edge in edges]
        if len(converted) < 2:
            raise RuntimeError(
                "[lepRecoSFProducer] %s has fewer than two bin edges." % label
            )
        for index in range(len(converted) - 1):
            if not converted[index + 1] > converted[index]:
                raise RuntimeError(
                    "[lepRecoSFProducer] Non-monotonic edges in %s: %r"
                    % (label, edges)
                )
        return converted

    @staticmethod
    def _find_category(node, key):
        if not isinstance(node, dict) or node.get("nodetype") != "category":
            raise RuntimeError(
                "[lepRecoSFProducer] Expected a category node, got %r"
                % (node.get("nodetype") if isinstance(node, dict) else type(node))
            )

        for item in node.get("content", []):
            if item.get("key") == key:
                return item.get("value")

        available = [item.get("key") for item in node.get("content", [])]
        raise RuntimeError(
            "[lepRecoSFProducer] Category key %r is unavailable for input %r. "
            "Available keys: %r"
            % (key, node.get("input"), available)
        )

    @staticmethod
    def _find_bin(edges, value, input_name, flow):
        value = float(value)
        index = bisect.bisect_right(edges, value) - 1

        if 0 <= index < len(edges) - 1:
            return index

        if flow == "clamp":
            if value < edges[0]:
                return 0
            return len(edges) - 2

        raise RuntimeError(
            "[lepRecoSFProducer] Value %s is outside the %s axis [%s, %s], "
            "with flow=%r."
            % (value, input_name, edges[0], edges[-1], flow)
        )

    def _load(self):
        if not os.path.exists(self.json_path):
            raise RuntimeError(
                "[lepRecoSFProducer] Electron reco SF JSON does not exist: %s"
                % self.json_path
            )

        payload = self._open_json(self.json_path)

        if int(payload.get("schema_version", -1)) != 2:
            raise RuntimeError(
                "[lepRecoSFProducer] Unsupported correction schema version %r in %s"
                % (payload.get("schema_version"), self.json_path)
            )

        matches = [
            correction
            for correction in payload.get("corrections", [])
            if correction.get("name") == self.correction_name
        ]
        if len(matches) != 1:
            raise RuntimeError(
                "[lepRecoSFProducer] Expected exactly one correction named %r "
                "in %s, found %d."
                % (self.correction_name, self.json_path, len(matches))
            )

        self.correction = matches[0]
        self.input_names = [
            item.get("name") for item in self.correction.get("inputs", [])
        ]

        expected = ["year", "ValType", "WorkingPoint", "eta", "pt"]
        if self.input_names != expected:
            raise RuntimeError(
                "[lepRecoSFProducer] Unexpected input schema for %s: %r; "
                "expected %r."
                % (self.correction_name, self.input_names, expected)
            )

        # Validate every multibinning node once at startup.
        self._validate_node(self.correction.get("data"), path=())

    def _validate_node(self, node, path):
        if not isinstance(node, dict):
            if not isinstance(node, (int, float)):
                raise RuntimeError(
                    "[lepRecoSFProducer] Unsupported leaf at %s: %r"
                    % ("/".join(path), type(node))
                )
            return

        nodetype = node.get("nodetype")

        if nodetype == "category":
            content = node.get("content", [])
            if not content:
                raise RuntimeError(
                    "[lepRecoSFProducer] Empty category node at %s"
                    % "/".join(path)
                )
            for item in content:
                self._validate_node(
                    item.get("value"),
                    path + ("%s=%s" % (node.get("input"), item.get("key")),),
                )
            return

        if nodetype == "multibinning":
            inputs = node.get("inputs", [])
            edge_sets = node.get("edges", [])
            if len(inputs) != len(edge_sets):
                raise RuntimeError(
                    "[lepRecoSFProducer] inputs/edges size mismatch at %s"
                    % "/".join(path)
                )

            n_total = 1
            for input_name, edges in zip(inputs, edge_sets):
                converted = self._validate_edges(
                    edges, "%s/%s" % ("/".join(path), input_name)
                )
                n_total *= len(converted) - 1

            content = node.get("content", [])
            if len(content) != n_total:
                raise RuntimeError(
                    "[lepRecoSFProducer] Multibinning content size mismatch at %s: "
                    "expected %d, got %d."
                    % ("/".join(path), n_total, len(content))
                )

            for value in content:
                if not isinstance(value, (int, float)):
                    raise RuntimeError(
                        "[lepRecoSFProducer] Non-numeric multibinning content at %s: %r"
                        % ("/".join(path), type(value))
                    )
            return

        raise RuntimeError(
            "[lepRecoSFProducer] Unsupported node type %r at %s"
            % (nodetype, "/".join(path))
        )

    def evaluate(self, year, valtype, working_point, eta, pt):
        node = self.correction["data"]
        node = self._find_category(node, str(year))
        node = self._find_category(node, str(valtype))
        node = self._find_category(node, str(working_point))

        if not isinstance(node, dict) or node.get("nodetype") != "multibinning":
            raise RuntimeError(
                "[lepRecoSFProducer] Expected multibinning after "
                "year/ValType/WorkingPoint."
            )

        inputs = node["inputs"]
        if inputs != ["eta", "pt"]:
            raise RuntimeError(
                "[lepRecoSFProducer] Unexpected multibinning inputs %r; "
                "expected ['eta', 'pt']." % inputs
            )

        eta_edges = self._validate_edges(node["edges"][0], "eta")
        pt_edges = self._validate_edges(node["edges"][1], "pt")
        flow = node.get("flow", "error")

        eta_bin = self._find_bin(eta_edges, eta, "eta", flow)
        pt_bin = self._find_bin(pt_edges, pt, "pt", flow)
        n_pt_bins = len(pt_edges) - 1
        flat_index = eta_bin * n_pt_bins + pt_bin

        value = float(node["content"][flat_index])
        if not math.isfinite(value) or value <= 0.0:
            raise RuntimeError(
                "[lepRecoSFProducer] Invalid value in %s: year=%s, "
                "ValType=%s, WP=%s, eta=%s, pt=%s, bins=(%d,%d), value=%s"
                % (
                    self.json_path,
                    year,
                    valtype,
                    working_point,
                    eta,
                    pt,
                    eta_bin,
                    pt_bin,
                    value,
                )
            )
        return value


class ElectronRecoSFProducer(Module):
    """
    Official EGM electron reconstruction SF producer for 2024 NanoAODv15.

    Applied working points:
      RecoBelow20 : 10 <= pT < 20 GeV
      Reco20to75  : 20 <= pT < 75 GeV
      RecoAbove75 : pT >= 75 GeV

    The analysis-specific promptMVA/selection SF is separate and is applied by
    lepSFProducer_2024_signedEta.py.
    """

    def __init__(self, *args, **kwargs):
        cmssw_base = os.environ.get("CMSSW_BASE", "")
        default_local = os.path.join(
            cmssw_base,
            "src", "TZWi", "TopAnalysis", "data", "leptonSF", "2024",
            "electron.json.gz",
        )
        default_cvmfs = (
            "/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/"
            "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/"
            "latest/electron.json.gz"
        )

        self.year = kwargs.get("year", "2024Prompt")
        self.correction_name = kwargs.get(
            "correctionName", "Electron-ID-SF"
        )
        self.max_match_dr = float(kwargs.get("maxMatchDR", 0.02))
        self.max_match_relpt = float(kwargs.get("maxMatchRelPt", 0.02))
        self.max_fallback_warnings = int(
            kwargs.get("maxFallbackWarnings", 5)
        )
        self._fallback_warnings = 0

        requested_json = kwargs.get("jsonPath")
        candidates = (
            [requested_json]
            if requested_json
            else [default_local, default_cvmfs]
        )

        self.evaluator = None
        self.json_path = None
        errors = []

        for candidate in candidates:
            if not candidate:
                continue
            try:
                print(
                    "[lepRecoSFProducer] Trying official 2024 electron reco "
                    "SF JSON with the built-in schema-v2 evaluator: %s"
                    % candidate
                )
                evaluator = _ElectronRecoJSONEvaluator(
                    candidate,
                    correction_name=self.correction_name,
                )
                self.evaluator = evaluator
                self.json_path = candidate
                break
            except Exception as exc:
                errors.append((candidate, str(exc)))

        if self.evaluator is None:
            lines = [
                "[lepRecoSFProducer] Failed to load the official EGM electron "
                "reconstruction SF JSON.",
                "[lepRecoSFProducer] Tried:",
            ]
            for path, error in errors:
                lines.append(
                    "[lepRecoSFProducer]   - %s : %s" % (path, error)
                )
            raise RuntimeError("\n".join(lines))

        print(
            "[lepRecoSFProducer] Loaded %s from %s"
            % (self.correction_name, self.json_path)
        )
        self._validate_evaluator()

    def _validate_evaluator(self):
        validation_points = (
            ("RecoBelow20", 0.1, 15.0),
            ("Reco20to75", 0.1, 30.0),
            ("RecoAbove75", 0.1, 100.0),
        )

        for working_point, eta, pt in validation_points:
            value = self.evaluator.evaluate(
                self.year, "sf", working_point, eta, pt
            )
            print(
                "[lepRecoSFProducer] Validation: %s, eta=%.3f, pt=%.1f "
                "-> SF=%.8f"
                % (working_point, eta, pt, value)
            )

    def beginJob(self):
        pass

    def endJob(self):
        self.evaluator = None

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        for prefix in ("Lepton_RecoSF", "Electron_RecoSF"):
            self.out.branch(prefix, "F")
            self.out.branch(prefix + "Up", "F")
            self.out.branch(prefix + "Down", "F")
            self.out.branch(prefix + "err", "F")

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
                "[lepRecoSFProducer] WARNING: no Electron collection match "
                "for selected electron pdgId=%d, pt=%.6g, eta=%.6g, "
                "phi=%.6g; using Lepton_eta as the reco-SF eta."
                % (int(pdgid), float(pt), float(eta), float(phi))
            )
            self._fallback_warnings += 1

        return fallback_eta

    @staticmethod
    def _working_point(pt):
        pt = float(pt)
        if pt < 20.0:
            return "RecoBelow20"
        if pt < 75.0:
            return "Reco20to75"
        return "RecoAbove75"

    def _evaluate(self, valtype, working_point, sc_eta, pt):
        # The official correction starts at 10 GeV.  Protect tiny floating
        # point excursions at the analysis threshold.
        pt_eval = max(float(pt), 10.000001)
        return self.evaluator.evaluate(
            self.year,
            valtype,
            working_point,
            float(sc_eta),
            pt_eval,
        )

    def _one_electron_sf(self, event, pdgid, pt, eta, phi):
        if abs(int(pdgid)) != 11:
            return 1.0, 1.0, 1.0

        sc_eta = self._matched_electron_sc_eta(
            event, pdgid, pt, eta, phi
        )
        working_point = self._working_point(pt)

        central = self._evaluate(
            "sf", working_point, sc_eta, pt
        )
        up = self._evaluate(
            "sfup", working_point, sc_eta, pt
        )
        down = self._evaluate(
            "sfdown", working_point, sc_eta, pt
        )

        return central, up, down

    @staticmethod
    def _product(values):
        result = 1.0
        for value in values:
            result *= float(value)
        return float(result)

    @staticmethod
    def _read_selected_lepton(event, index):
        tree = event._tree
        return (
            int(getattr(tree, "b_out_Lepton%d_pdgId" % index)),
            float(getattr(tree, "b_out_Lepton%d_pt" % index)),
            float(getattr(tree, "b_out_Lepton%d_eta" % index)),
            float(getattr(tree, "b_out_Lepton%d_phi" % index)),
        )

    def analyze(self, event):
        central_values = []
        up_values = []
        down_values = []

        for index in (1, 2, 3):
            pdgid, pt, eta, phi = self._read_selected_lepton(
                event, index
            )
            central, up, down = self._one_electron_sf(
                event, pdgid, pt, eta, phi
            )
            central_values.append(central)
            up_values.append(up)
            down_values.append(down)

        central_total = self._product(central_values)
        up_total = self._product(up_values)
        down_total = self._product(down_values)
        error_total = max(
            abs(up_total - central_total),
            abs(central_total - down_total),
        )

        for prefix in ("Lepton_RecoSF", "Electron_RecoSF"):
            self.out.fillBranch(prefix, central_total)
            self.out.fillBranch(prefix + "Up", up_total)
            self.out.fillBranch(prefix + "Down", down_total)
            self.out.fillBranch(prefix + "err", error_total)

        return True


lepRecoSF = lambda: ElectronRecoSFProducer()
lepRecoSF_2024 = lambda: ElectronRecoSFProducer()
