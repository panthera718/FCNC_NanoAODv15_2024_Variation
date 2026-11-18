import os
from math import sqrt, pi

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class _HistogramWeightCalculator(object):
    """
    Read a value or its bin error from a TH1/TH2.

    Underflow and overflow are clamped to the nearest valid edge bin.
    For the 2024 promptMVA lepton SF maps used here, the TH2 axes are:
      x-axis = signed eta or signed SC eta
      y-axis = pT
    """

    def __init__(self, histogram):
        self.histogram = histogram

    def _clamp_xbin(self, x):
        axis = self.histogram.GetXaxis()
        return max(1, min(self.histogram.GetNbinsX(), axis.FindBin(float(x))))

    def _clamp_ybin(self, y):
        axis = self.histogram.GetYaxis()
        return max(1, min(self.histogram.GetNbinsY(), axis.FindBin(float(y))))

    def getWeight(self, x, y=0.0):
        if self.histogram is None:
            raise RuntimeError("The input SF histogram is not loaded.")
        if not self.histogram.InheritsFrom("TH2"):
            return float(self.histogram.GetBinContent(self._clamp_xbin(x)))
        return float(self.histogram.GetBinContent(self._clamp_xbin(x), self._clamp_ybin(y)))

    def getWeightErr(self, x, y=0.0):
        if self.histogram is None:
            raise RuntimeError("The input SF histogram is not loaded.")
        if not self.histogram.InheritsFrom("TH2"):
            return float(self.histogram.GetBinError(self._clamp_xbin(x)))
        return float(self.histogram.GetBinError(self._clamp_xbin(x), self._clamp_ybin(y)))


class _PromptMVALeptonSFMap(object):
    """
    Small pure-Python SF reader for the analysis-specific Run-3 promptMVA SF maps.

    The maps produced by make_lepton_sf_root_lumiNorm_skipBad_signedEta.py store
    the central Data/MC SF in h2_sf.  The current map contains the central value
    and its bin error only; no separate systematic histograms are assumed here.
    """

    def __init__(self, root_file_name, hist_name="h2_sf", label=""):
        self.root_file_name = root_file_name
        self.hist_name = hist_name
        self.label = label if label else hist_name
        self.hist = None
        self._load()

    def _load(self):
        root_file = ROOT.TFile.Open(self.root_file_name, "READ")
        if not root_file or root_file.IsZombie():
            raise RuntimeError("Cannot open lepton SF ROOT file: {0}".format(self.root_file_name))

        obj = root_file.Get(self.hist_name)
        if not obj:
            root_file.Close()
            raise RuntimeError(
                "Histogram {0} is not found in lepton SF ROOT file: {1}".format(
                    self.hist_name, self.root_file_name
                )
            )

        self.hist = obj.Clone("lepton_sf_{0}_{1}".format(self.label, self.hist_name).replace("/", "_"))
        self.hist.SetDirectory(0)
        root_file.Close()

        if not self.hist.InheritsFrom("TH2"):
            raise RuntimeError(
                "Histogram {0} in {1} is not a TH2. This producer expects a 2D eta-pT SF map.".format(
                    self.hist_name, self.root_file_name
                )
            )

        print("[lepSFProducer] Loaded {0} from {1}".format(self.hist_name, self.root_file_name))
        print("[lepSFProducer]   x-axis: {0} to {1}, nbins={2}".format(
            self.hist.GetXaxis().GetXmin(), self.hist.GetXaxis().GetXmax(), self.hist.GetNbinsX()
        ))
        print("[lepSFProducer]   y-axis: {0} to {1}, nbins={2}".format(
            self.hist.GetYaxis().GetXmin(), self.hist.GetYaxis().GetXmax(), self.hist.GetNbinsY()
        ))

    def getSF(self, eta, pt):
        wc = _HistogramWeightCalculator(self.hist)
        return wc.getWeight(float(eta), float(pt))

    def getSFErr(self, eta, pt):
        wc = _HistogramWeightCalculator(self.hist)
        return wc.getWeightErr(float(eta), float(pt))


class lepSFProducer(Module):
    """
    Apply the analysis-specific 2024 Run-3 promptMVA lepton SFs.

    Input:
      - This module must run after FCNCTriLepton, because it reads
        b_out_Lepton{1,2,3}_{pdgId,pt,eta,phi} cached by that module.
      - The SF maps are the signed-eta maps:
          data/leptonSF/leptonSF_electron_2024_signedEta.root
          data/leptonSF/leptonSF_muon_2024_signedEta.root

    Eta convention:
      - Electron lookup uses signed supercluster eta when the original NanoAOD
        Electron object can be matched:
            scEta = Electron_eta + Electron_deltaEtaSC
        If matching fails, it falls back to the signed LeptonN_eta from FCNCTriLepton.
      - Muon lookup uses signed Muon eta.
      - No abs(eta) is applied.
    """

    def __init__(
        self,
        electron_file=None,
        muon_file=None,
        hist_name="h2_sf",
    ):
        base = os.path.join(os.environ["CMSSW_BASE"], "src", "TZWi", "TopAnalysis", "data", "leptonSF")

        if electron_file is None:
            electron_file = os.path.join(base, "leptonSF_electron_2024_signedEta.root")
        if muon_file is None:
            muon_file = os.path.join(base, "leptonSF_muon_2024_signedEta.root")

        self.electron_file = electron_file
        self.muon_file = muon_file
        self.hist_name = hist_name

        self._worker_el = None
        self._worker_mu = None

    def beginJob(self):
        self._worker_el = _PromptMVALeptonSFMap(self.electron_file, self.hist_name, "electron")
        self._worker_mu = _PromptMVALeptonSFMap(self.muon_file, self.hist_name, "muon")

    def endJob(self):
        if hasattr(self, "_worker_el"):
            del self._worker_el
        if hasattr(self, "_worker_mu"):
            del self._worker_mu

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Lepton_SF", "F")
        self.out.branch("Lepton_SFerr", "F")
        self.out.branch("Lepton_SFerrStat", "F")
        self.out.branch("Lepton_SFerrSyst", "F")
        self.out.branch("Electron_SF", "F")
        self.out.branch("Electron_SFerr", "F")
        self.out.branch("Electron_SFerrStat", "F")
        self.out.branch("Electron_SFerrSyst", "F")
        self.out.branch("Muon_SF", "F")
        self.out.branch("Muon_SFerr", "F")
        self.out.branch("Muon_SFerrStat", "F")
        self.out.branch("Muon_SFerrSyst", "F")

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

    def _matched_electron_sc_eta(self, event, pdgid, pt, eta, phi):
        """
        Return signed electron SC eta by matching the selected LeptonN back to
        the original NanoAOD Electron collection.

        FCNCTriLepton does not write LeptonN_scEta.  Therefore, for electrons,
        use the original Electron_deltaEtaSC if a unique tight numerical match
        can be found.  If not, return the selected lepton eta as a safe fallback.
        """
        fallback_eta = float(eta)

        try:
            electrons = Collection(event, "Electron")
        except Exception:
            return fallback_eta

        target_pdgid = int(pdgid)
        best = None
        best_score = None

        for ele in electrons:
            try:
                ele_pdgid = -11 * int(ele.charge)
                if ele_pdgid != target_pdgid:
                    continue

                dr2 = self._delta_r2(eta, phi, ele.eta, ele.phi)
                dpt_rel = abs(float(pt) - float(ele.pt)) / max(1.0, abs(float(pt)))
                score = dr2 + dpt_rel * dpt_rel

                if best_score is None or score < best_score:
                    best_score = score
                    best = ele
            except Exception:
                continue

        # The selected lepton p4 is built directly from the NanoAOD electron p4,
        # so a genuine match should be extremely close.  Keep a loose fallback
        # guard to avoid accidentally matching the wrong electron in rare cases.
        if best is None or best_score is None or best_score > 1.0e-6:
            return fallback_eta

        try:
            return float(best.eta) + float(best.deltaEtaSC)
        except Exception:
            return float(best.eta)

    @staticmethod
    def _combine_product_with_errors(components):
        """
        components: list of (sf, total_error, stat_error, syst_error)

        For a product of independent SFs, propagate uncertainties through
        relative errors:
            sigma(product) = product * sqrt(sum_i (sigma_i / sf_i)^2)
        """
        product = 1.0
        for sf, _, _, _ in components:
            product *= float(sf)

        def combine(which):
            rel2 = 0.0
            abs2_if_zero = 0.0
            for sf, err_total, err_stat, err_syst in components:
                if which == "total":
                    err = float(err_total)
                elif which == "stat":
                    err = float(err_stat)
                elif which == "syst":
                    err = float(err_syst)
                else:
                    err = 0.0

                sf = float(sf)
                if sf != 0.0:
                    rel2 += (err / sf) ** 2.0
                else:
                    abs2_if_zero += err * err

            if product != 0.0:
                return abs(product) * sqrt(max(0.0, rel2))
            return sqrt(max(0.0, abs2_if_zero))

        return product, combine("total"), combine("stat"), combine("syst")

    def _process_one_lepton(self, event, idx):
        pdgid = int(getattr(event._tree, "b_out_Lepton%d_pdgId" % idx))
        pt = float(getattr(event._tree, "b_out_Lepton%d_pt" % idx))
        eta = float(getattr(event._tree, "b_out_Lepton%d_eta" % idx))
        phi = float(getattr(event._tree, "b_out_Lepton%d_phi" % idx))

        apdgid = abs(pdgid)
        if apdgid == 11:
            lookup_eta = self._matched_electron_sc_eta(event, pdgid, pt, eta, phi)
            sf = self._worker_el.getSF(lookup_eta, pt)
            sf_err = self._worker_el.getSFErr(lookup_eta, pt)

            # The current h2_sf bin error is treated as the statistical error.
            # No dedicated systematic uncertainty map is available yet.
            return "electron", sf, sf_err, sf_err, 0.0

        if apdgid == 13:
            # The muon SF map is signed-eta, so do not use abs(eta).
            lookup_eta = eta
            sf = self._worker_mu.getSF(lookup_eta, pt)
            sf_err = self._worker_mu.getSFErr(lookup_eta, pt)
            return "muon", sf, sf_err, sf_err, 0.0

        return None, 1.0, 0.0, 0.0, 0.0

    def analyze(self, event):
        all_components = []
        electron_components = []
        muon_components = []

        for idx in (1, 2, 3):
            flavor, sf, sf_err, sf_err_stat, sf_err_syst = self._process_one_lepton(event, idx)
            component = (sf, sf_err, sf_err_stat, sf_err_syst)

            all_components.append(component)
            if flavor == "electron":
                electron_components.append(component)
            elif flavor == "muon":
                muon_components.append(component)

        sf_lepton, sf_lepton_err, sf_lepton_err_stat, sf_lepton_err_syst = self._combine_product_with_errors(all_components)
        sf_el, sf_el_err, sf_el_err_stat, sf_el_err_syst = self._combine_product_with_errors(electron_components)
        sf_mu, sf_mu_err, sf_mu_err_stat, sf_mu_err_syst = self._combine_product_with_errors(muon_components)

        self.out.fillBranch("Lepton_SF", sf_lepton)
        self.out.fillBranch("Lepton_SFerr", sf_lepton_err)
        self.out.fillBranch("Lepton_SFerrStat", sf_lepton_err_stat)
        self.out.fillBranch("Lepton_SFerrSyst", sf_lepton_err_syst)

        self.out.fillBranch("Electron_SF", sf_el)
        self.out.fillBranch("Electron_SFerr", sf_el_err)
        self.out.fillBranch("Electron_SFerrStat", sf_el_err_stat)
        self.out.fillBranch("Electron_SFerrSyst", sf_el_err_syst)

        self.out.fillBranch("Muon_SF", sf_mu)
        self.out.fillBranch("Muon_SFerr", sf_mu_err)
        self.out.fillBranch("Muon_SFerrStat", sf_mu_err_stat)
        self.out.fillBranch("Muon_SFerrSyst", sf_mu_err_syst)

        return True


# Define modules using the syntax "name = lambda : constructor" to avoid
# loading the ROOT files unless the module is explicitly requested.
lepSF = lambda: lepSFProducer()
