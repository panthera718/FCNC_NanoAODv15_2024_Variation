import os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class TauTriggerSFProducer(Module, object):
    def __init__(self, *args, **kwargs):
        self.mode = kwargs.get("mode")
        self.year = kwargs.get("year", "2018")

        # normal modes use tight tau WP, NPL modes use fakeable tau WP
        self.wp_tight = kwargs.get("wp_tight", "Medium")
        self.wp_fakeable = kwargs.get("wp_fakeable", "VVLoose")

        if self.year != "2018":
            raise RuntimeError("This producer is currently configured only for year=2018")

        self.root_path = os.path.join(
            os.environ["CMSSW_BASE"],
            "src",
            "TZWi",
            "TopAnalysis",
            "data",
            "tauTriggerSF",
            "2018UL_tauTriggerEff_DeepTau2017v2p1_NanoAOD.root",
        )

        self.tf = None
        self._hist_cache = {}

    def beginJob(self):
        self.tf = ROOT.TFile.Open(self.root_path)
        if not self.tf or self.tf.IsZombie():
            raise RuntimeError("Cannot open tau trigger SF ROOT file: %s" % self.root_path)

    def endJob(self):
        if self.tf:
            self.tf.Close()
            self.tf = None

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch("TauTriggerWeight", "F")
        self.out.branch("TauTriggerWeight_up", "F")
        self.out.branch("TauTriggerWeight_down", "F")

        self.out.branch("TauTriggerEffMC", "F")
        self.out.branch("TauTriggerEffData", "F")
        self.out.branch("TauTriggerEffData_up", "F")
        self.out.branch("TauTriggerEffData_down", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    @staticmethod
    def _set_tree_attr(tree, name, value):
        setattr(tree, "b_out_%s" % name, value)

    def _base_mode(self):
        if self.mode is None:
            return ""
        if self.mode.startswith("NPL"):
            return self.mode[3:]
        return self.mode

    def _is_npl_mode(self):
        return self.mode is not None and self.mode.startswith("NPL")

    def _active_wp(self):
        if self._is_npl_mode():
            return self.wp_fakeable
        return self.wp_tight

    def _get_hist(self, name):
        if name in self._hist_cache:
            return self._hist_cache[name]

        h = self.tf.Get(name)
        if not h:
            raise RuntimeError("Histogram not found: %s in %s" % (name, self.root_path))

        h.SetDirectory(0)
        self._hist_cache[name] = h
        return h

    @staticmethod
    def _clip01(x):
        if x < 0.0:
            return 0.0
        if x > 1.0:
            return 1.0
        return x

    def _get_bin_and_content(self, hist, x):
        xaxis = hist.GetXaxis()
        xx = float(x)

        xmin = xaxis.GetXmin()
        xmax = xaxis.GetXmax()

        if xx < xmin:
            xx = xmin + 1e-6
        if xx >= xmax:
            xx = xmax - 1e-6

        ibin = xaxis.FindFixBin(xx)
        if ibin < 1:
            ibin = 1
        if ibin > hist.GetNbinsX():
            ibin = hist.GetNbinsX()

        return ibin, float(hist.GetBinContent(ibin)), float(hist.GetBinError(ibin))

    @staticmethod
    def _dm_label(dm):
        dm = int(dm)
        if dm == 0:
            return "dm0"
        if dm == 1:
            return "dm1"
        if dm == 10:
            return "dm10"
        if dm == 11:
            return "dm11"
        return "dmall"

    def _hist_name(self, kind, trig, tau):
        # kind = "mc" or "sf"
        # trig = "etau", "mutau", "ditau"
        dm_label = self._dm_label(tau["dm"])
        wp = self._active_wp()
        return "%s_%s_%s_%s_fitted" % (kind, trig, wp, dm_label)

    def _tau_leg_mc_eff(self, trig, tau):
        # central tau-trigger SF/eff is applied only to genuine hadronic tau
        if int(tau["gen"]) != 5:
            return 0.0

        h = self._get_hist(self._hist_name("mc", trig, tau))
        _, eff, _ = self._get_bin_and_content(h, tau["pt"])
        return self._clip01(eff)

    def _tau_leg_sf_and_err(self, trig, tau):
        if int(tau["gen"]) != 5:
            return 1.0, 0.0

        h = self._get_hist(self._hist_name("sf", trig, tau))
        _, sf, err = self._get_bin_and_content(h, tau["pt"])

        if sf < 0.0:
            sf = 0.0
        if err < 0.0:
            err = 0.0

        return sf, err

    def _tau_leg_data_eff(self, trig, tau, variation="central"):
        """
        variation = "central", "up", "down"
        data_eff = mc_eff * sf_var
        """
        mc_eff = self._tau_leg_mc_eff(trig, tau)
        sf, err = self._tau_leg_sf_and_err(trig, tau)

        if variation == "up":
            sf_var = sf + err
        elif variation == "down":
            sf_var = max(0.0, sf - err)
        else:
            sf_var = sf

        return self._clip01(mc_eff * sf_var)

    @staticmethod
    def _prob_or(prob_list):
        p_not = 1.0
        for p in prob_list:
            pp = p
            if pp < 0.0:
                pp = 0.0
            if pp > 1.0:
                pp = 1.0
            p_not *= (1.0 - pp)
        return 1.0 - p_not

    def _single_tau_cross_eff(self, taus, trig, kind="mc"):
        """
        For etau or mutau tau-leg event efficiency:
          event eff = OR over all selected taus
        kind = "mc", "data", "data_up", "data_down"
        """
        leg_effs = []

        for tau in taus:
            if kind == "mc":
                leg_effs.append(self._tau_leg_mc_eff(trig, tau))
            elif kind == "data":
                leg_effs.append(self._tau_leg_data_eff(trig, tau, "central"))
            elif kind == "data_up":
                leg_effs.append(self._tau_leg_data_eff(trig, tau, "up"))
            elif kind == "data_down":
                leg_effs.append(self._tau_leg_data_eff(trig, tau, "down"))
            else:
                raise RuntimeError("Unknown kind: %s" % kind)

        return self._prob_or(leg_effs)

    def _ditau_eff(self, taus, kind="mc"):
        """
        Factorized tau-leg approximation for ditau:
          - if 2 taus: e1*e2
          - if 3 taus: OR over all tau pairs
        kind = "mc", "data", "data_up", "data_down"
        """
        n = len(taus)
        if n < 2:
            return 0.0

        leg_effs = []
        for tau in taus:
            if kind == "mc":
                leg_effs.append(self._tau_leg_mc_eff("ditau", tau))
            elif kind == "data":
                leg_effs.append(self._tau_leg_data_eff("ditau", tau, "central"))
            elif kind == "data_up":
                leg_effs.append(self._tau_leg_data_eff("ditau", tau, "up"))
            elif kind == "data_down":
                leg_effs.append(self._tau_leg_data_eff("ditau", tau, "down"))
            else:
                raise RuntimeError("Unknown kind: %s" % kind)

        pair_probs = []
        for i in range(n):
            for j in range(i + 1, n):
                pair_probs.append(leg_effs[i] * leg_effs[j])

        return self._prob_or(pair_probs)

    @staticmethod
    def _combine_or(a, b):
        aa = a
        bb = b
        if aa < 0.0:
            aa = 0.0
        if aa > 1.0:
            aa = 1.0
        if bb < 0.0:
            bb = 0.0
        if bb > 1.0:
            bb = 1.0
        return 1.0 - (1.0 - aa) * (1.0 - bb)

    def _get_selected_taus(self, tree):
        taus = []
        for idx in (1, 2, 3):
            pdgid_name = "b_out_Lepton%d_pdgId" % idx
            if not hasattr(tree, pdgid_name):
                continue

            pdgid = abs(int(getattr(tree, pdgid_name)))
            if pdgid != 15:
                continue

            pt = float(getattr(tree, "b_out_Lepton%d_pt" % idx))
            dm = int(getattr(tree, "b_out_Lepton%d_tauDecayMode" % idx))

            gen_name = "b_out_Lepton%d_tauGenPartFlav" % idx
            gen = int(getattr(tree, gen_name)) if hasattr(tree, gen_name) else -1

            taus.append({
                "pt": pt,
                "dm": dm,
                "gen": gen,
            })

        return taus

    def _compute_efficiencies(self, tree):
        """
        Returns:
          mc_eff_central, data_eff_central, data_eff_up, data_eff_down
        """
        mode = self._base_mode()
        taus = self._get_selected_taus(tree)

        # default for channels with no tau-trigger application
        mc_eff = 1.0
        data_eff = 1.0
        data_eff_up = 1.0
        data_eff_down = 1.0

        if mode == "TaElEl":
            mc_eff = self._single_tau_cross_eff(taus, "etau", kind="mc")
            data_eff = self._single_tau_cross_eff(taus, "etau", kind="data")
            data_eff_up = self._single_tau_cross_eff(taus, "etau", kind="data_up")
            data_eff_down = self._single_tau_cross_eff(taus, "etau", kind="data_down")

        elif mode == "TaMuMu":
            mc_eff = self._single_tau_cross_eff(taus, "mutau", kind="mc")
            data_eff = self._single_tau_cross_eff(taus, "mutau", kind="data")
            data_eff_up = self._single_tau_cross_eff(taus, "mutau", kind="data_up")
            data_eff_down = self._single_tau_cross_eff(taus, "mutau", kind="data_down")

        elif mode == "ElTaTa":
            mc_etau = self._single_tau_cross_eff(taus, "etau", kind="mc")
            data_etau = self._single_tau_cross_eff(taus, "etau", kind="data")
            data_etau_up = self._single_tau_cross_eff(taus, "etau", kind="data_up")
            data_etau_down = self._single_tau_cross_eff(taus, "etau", kind="data_down")

            mc_ditau = self._ditau_eff(taus, kind="mc")
            data_ditau = self._ditau_eff(taus, kind="data")
            data_ditau_up = self._ditau_eff(taus, kind="data_up")
            data_ditau_down = self._ditau_eff(taus, kind="data_down")

            mc_eff = self._combine_or(mc_etau, mc_ditau)
            data_eff = self._combine_or(data_etau, data_ditau)
            data_eff_up = self._combine_or(data_etau_up, data_ditau_up)
            data_eff_down = self._combine_or(data_etau_down, data_ditau_down)

        elif mode == "MuTaTa":
            mc_mutau = self._single_tau_cross_eff(taus, "mutau", kind="mc")
            data_mutau = self._single_tau_cross_eff(taus, "mutau", kind="data")
            data_mutau_up = self._single_tau_cross_eff(taus, "mutau", kind="data_up")
            data_mutau_down = self._single_tau_cross_eff(taus, "mutau", kind="data_down")

            mc_ditau = self._ditau_eff(taus, kind="mc")
            data_ditau = self._ditau_eff(taus, kind="data")
            data_ditau_up = self._ditau_eff(taus, kind="data_up")
            data_ditau_down = self._ditau_eff(taus, kind="data_down")

            mc_eff = self._combine_or(mc_mutau, mc_ditau)
            data_eff = self._combine_or(data_mutau, data_ditau)
            data_eff_up = self._combine_or(data_mutau_up, data_ditau_up)
            data_eff_down = self._combine_or(data_mutau_down, data_ditau_down)

        elif mode == "TaTaTa":
            mc_eff = self._ditau_eff(taus, kind="mc")
            data_eff = self._ditau_eff(taus, kind="data")
            data_eff_up = self._ditau_eff(taus, kind="data_up")
            data_eff_down = self._ditau_eff(taus, kind="data_down")

        return float(mc_eff), float(data_eff), float(data_eff_up), float(data_eff_down)

    def analyze(self, event):
        tree = event._tree

        mc_eff, data_eff, data_eff_up, data_eff_down = self._compute_efficiencies(tree)

        if mc_eff <= 0.0:
            weight = 1.0
            weight_up = 1.0
            weight_down = 1.0
        else:
            weight = data_eff / mc_eff
            weight_up = data_eff_up / mc_eff
            weight_down = data_eff_down / mc_eff

        self._set_tree_attr(tree, "TauTriggerWeight", float(weight))
        self._set_tree_attr(tree, "TauTriggerWeight_up", float(weight_up))
        self._set_tree_attr(tree, "TauTriggerWeight_down", float(weight_down))

        self._set_tree_attr(tree, "TauTriggerEffMC", float(mc_eff))
        self._set_tree_attr(tree, "TauTriggerEffData", float(data_eff))
        self._set_tree_attr(tree, "TauTriggerEffData_up", float(data_eff_up))
        self._set_tree_attr(tree, "TauTriggerEffData_down", float(data_eff_down))

        self.out.fillBranch("TauTriggerWeight", float(weight))
        self.out.fillBranch("TauTriggerWeight_up", float(weight_up))
        self.out.fillBranch("TauTriggerWeight_down", float(weight_down))

        self.out.fillBranch("TauTriggerEffMC", float(mc_eff))
        self.out.fillBranch("TauTriggerEffData", float(data_eff))
        self.out.fillBranch("TauTriggerEffData_up", float(data_eff_up))
        self.out.fillBranch("TauTriggerEffData_down", float(data_eff_down))

        return True


def _factory(mode, year="2018", wp_tight="Medium", wp_fakeable="VVLoose"):
    return lambda: TauTriggerSFProducer(
        mode=mode,
        year=year,
        wp_tight=wp_tight,
        wp_fakeable=wp_fakeable,
    )


MODES = [
    "ElElEl", "MuElEl", "ElMuMu", "MuMuMu",
    "TaTaTa", "ElTaTa", "MuTaTa", "TaElEl", "TaMuMu",
    "NPLElElEl", "NPLMuElEl", "NPLElMuMu", "NPLMuMuMu",
    "NPLTaTaTa", "NPLElTaTa", "NPLMuTaTa", "NPLTaElEl", "NPLTaMuMu",
]

for _mode in MODES:
    globals()["tauTriggerSF_%s_2018" % _mode] = _factory(
        _mode,
        year="2018",
        wp_tight="Medium",
        wp_fakeable="VVLoose",
    )