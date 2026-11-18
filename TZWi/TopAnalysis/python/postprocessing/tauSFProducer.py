import os
import math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from correctionlib import _core
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class TauSFProducer(Module, object):
    def __init__(self, *args, **kwargs):
        self.mode = kwargs.get("mode")
        self.year = kwargs.get("year", "2018")

        # Exact supported nominal WP combination in the provided tau.json
        self.vsjet_wp = kwargs.get("vsjet_wp", "Medium")
        self.vsjet_wp_vse = kwargs.get("vsjet_wp_vse", "VVLoose")
        self.vsjet_flag = kwargs.get("vsjet_flag", "pt")

        self.vse_wp = kwargs.get("vse_wp", "VVLoose")
        self.vsmu_wp = kwargs.get("vsmu_wp", "VLoose")

        # Additional uncertainty recommended when using looser anti-lepton WPs
        # than the standard VSjet SF measurement setup.
        self.apply_other_vslepwp_unc = bool(kwargs.get("apply_other_vslepwp_unc", True))

        if self.year != "2018":
            raise RuntimeError("This producer is currently configured only for year=2018")

        self.json_path = os.path.join(
            os.environ["CMSSW_BASE"],
            "src",
            "TZWi",
            "TopAnalysis",
            "data",
            "tauSF",
            "tau.json.gz",
        )

        self.cset = None
        self.corr_vse = None
        self.corr_vsmu = None
        self.corr_vsjet = None

    def beginJob(self):
        self.cset = _core.CorrectionSet.from_file(self.json_path)

        # Exact correction names in tau.json
        self.corr_vse = self.cset["DeepTau2017v2p1VSe"]
        self.corr_vsmu = self.cset["DeepTau2017v2p1VSmu"]
        self.corr_vsjet = self.cset["DeepTau2017v2p1VSjet"]

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch("TauIDSF", "F")
        self.out.branch("TauIDSF_Lepton1", "F")
        self.out.branch("TauIDSF_Lepton2", "F")
        self.out.branch("TauIDSF_Lepton3", "F")

        self.out.branch("TauIDSF_up", "F")
        self.out.branch("TauIDSF_down", "F")

        self.out.branch("TauIDSF_genuineTau_up", "F")
        self.out.branch("TauIDSF_genuineTau_down", "F")

        self.out.branch("TauIDSF_eToTau_up", "F")
        self.out.branch("TauIDSF_eToTau_down", "F")

        self.out.branch("TauIDSF_muToTau_up", "F")
        self.out.branch("TauIDSF_muToTau_down", "F")

        self.out.branch("TauIDSF_otherVSlepWP_up", "F")
        self.out.branch("TauIDSF_otherVSlepWP_down", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    @staticmethod
    def _set_tree_attr(tree, name, value):
        setattr(tree, "b_out_%s" % name, value)

    @staticmethod
    def _prod(vals):
        out = 1.0
        for v in vals:
            out *= float(v)
        return float(out)

    @staticmethod
    def _is_tau(tree, idx):
        name = "b_out_Lepton%d_pdgId" % idx
        if not hasattr(tree, name):
            return False
        return abs(int(getattr(tree, name))) == 15

    @staticmethod
    def _tau_info(tree, idx):
        is_tau = TauSFProducer._is_tau(tree, idx)
        if not is_tau:
            return {
                "idx": idx,
                "is_tau": False,
                "pt": 0.0,
                "eta": 0.0,
                "dm": -1,
                "gen": -1,
            }

        return {
            "idx": idx,
            "is_tau": True,
            "pt": float(getattr(tree, "b_out_Lepton%d_pt" % idx)),
            "eta": float(getattr(tree, "b_out_Lepton%d_eta" % idx)),
            "dm": int(getattr(tree, "b_out_Lepton%d_tauDecayMode" % idx)),
            "gen": int(getattr(tree, "b_out_Lepton%d_tauGenPartFlav" % idx)),
        }

    @staticmethod
    def _other_vslepwp_rel_unc(pt):
        # Tau POG recommendation for looser VSe/VSmu than measurement setup:
        # 3% for pT < 100 GeV, 15% for pT > 100 GeV
        if float(pt) > 100.0:
            return 0.15
        return 0.03

    def _eval_vse(self, eta, genmatch, syst):
        # Exact signature in tau.json:
        #   DeepTau2017v2p1VSe.evaluate(eta, genmatch, wp, syst)
        return float(self.corr_vse.evaluate(float(eta), int(genmatch), self.vse_wp, syst))

    def _eval_vsmu(self, eta, genmatch, syst):
        # Exact signature in tau.json:
        #   DeepTau2017v2p1VSmu.evaluate(eta, genmatch, wp, syst)
        return float(self.corr_vsmu.evaluate(float(eta), int(genmatch), self.vsmu_wp, syst))

    def _eval_vsjet(self, pt, dm, genmatch, syst):
        # Exact signature in tau.json:
        #   DeepTau2017v2p1VSjet.evaluate(pt, dm, genmatch, wp, wp_VSe, syst, flag)
        #
        # In this tau.json, the nominal VSjet value is stored in the category "default"
        # rather than an explicit "nom" key. Passing "nom" here intentionally falls back
        # to that default nominal value.
        return float(
            self.corr_vsjet.evaluate(
                float(pt),
                int(dm),
                int(genmatch),
                self.vsjet_wp,
                self.vsjet_wp_vse,
                syst,
                self.vsjet_flag,
            )
        )

    def _sf_nominal_one(self, tau):
        if not tau["is_tau"]:
            return 1.0

        gen = int(tau["gen"])

        # genuine hadronic tau
        if gen == 5:
            return self._eval_vsjet(tau["pt"], tau["dm"], gen, "nom")

        # e -> tau fake
        if gen in (1, 3):
            return self._eval_vse(tau["eta"], gen, "nom")

        # mu -> tau fake
        if gen in (2, 4):
            return self._eval_vsmu(tau["eta"], gen, "nom")

        # jet fake / unmatched
        return 1.0

    def _event_weight_with_component(self, taus, component=None, direction=None):
        vals = []

        for tau in taus:
            nom = self._sf_nominal_one(tau)

            if component is None:
                vals.append(nom)
                continue

            gen = int(tau["gen"])

            if component == "genuineTau" and gen == 5:
                vals.append(self._eval_vsjet(tau["pt"], tau["dm"], gen, direction))
                continue

            if component == "eToTau" and gen in (1, 3):
                vals.append(self._eval_vse(tau["eta"], gen, direction))
                continue

            if component == "muToTau" and gen in (2, 4):
                vals.append(self._eval_vsmu(tau["eta"], gen, direction))
                continue

            # extra uncertainty recommended when using looser VSe/VSmu WPs
            # This applies to the genuine-tau VSjet SF.
            if component == "otherVSlepWP" and self.apply_other_vslepwp_unc and gen == 5:
                rel = self._other_vslepwp_rel_unc(tau["pt"])
                if direction == "up":
                    vals.append(nom * (1.0 + rel))
                elif direction == "down":
                    vals.append(max(0.0, nom * (1.0 - rel)))
                else:
                    vals.append(nom)
                continue

            vals.append(nom)

        return self._prod(vals)

    @staticmethod
    def _combine_total_from_components(nominal, component_pairs):
        if nominal <= 0.0:
            return 1.0, 1.0

        rel_up_sq = 0.0
        rel_dn_sq = 0.0

        for w_up, w_dn in component_pairs:
            rel_up = max(0.0, float(w_up) / nominal - 1.0)
            rel_dn = max(0.0, 1.0 - float(w_dn) / nominal)

            rel_up_sq += rel_up * rel_up
            rel_dn_sq += rel_dn * rel_dn

        total_up = nominal * (1.0 + math.sqrt(rel_up_sq))
        total_dn = max(0.0, nominal * (1.0 - math.sqrt(rel_dn_sq)))
        return float(total_up), float(total_dn)

    def analyze(self, event):
        tree = event._tree
        taus = [self._tau_info(tree, i) for i in (1, 2, 3)]

        sf1 = self._sf_nominal_one(taus[0])
        sf2 = self._sf_nominal_one(taus[1])
        sf3 = self._sf_nominal_one(taus[2])

        nominal = float(sf1 * sf2 * sf3)

        w_g_up = self._event_weight_with_component(taus, "genuineTau", "up")
        w_g_dn = self._event_weight_with_component(taus, "genuineTau", "down")

        w_e_up = self._event_weight_with_component(taus, "eToTau", "up")
        w_e_dn = self._event_weight_with_component(taus, "eToTau", "down")

        w_m_up = self._event_weight_with_component(taus, "muToTau", "up")
        w_m_dn = self._event_weight_with_component(taus, "muToTau", "down")

        w_o_up = self._event_weight_with_component(taus, "otherVSlepWP", "up")
        w_o_dn = self._event_weight_with_component(taus, "otherVSlepWP", "down")

        w_tot_up, w_tot_dn = self._combine_total_from_components(
            nominal,
            [
                (w_g_up, w_g_dn),
                (w_e_up, w_e_dn),
                (w_m_up, w_m_dn),
                (w_o_up, w_o_dn),
            ],
        )

        self._set_tree_attr(tree, "TauIDSF", nominal)
        self._set_tree_attr(tree, "TauIDSF_Lepton1", float(sf1))
        self._set_tree_attr(tree, "TauIDSF_Lepton2", float(sf2))
        self._set_tree_attr(tree, "TauIDSF_Lepton3", float(sf3))

        self._set_tree_attr(tree, "TauIDSF_up", w_tot_up)
        self._set_tree_attr(tree, "TauIDSF_down", w_tot_dn)

        self._set_tree_attr(tree, "TauIDSF_genuineTau_up", w_g_up)
        self._set_tree_attr(tree, "TauIDSF_genuineTau_down", w_g_dn)

        self._set_tree_attr(tree, "TauIDSF_eToTau_up", w_e_up)
        self._set_tree_attr(tree, "TauIDSF_eToTau_down", w_e_dn)

        self._set_tree_attr(tree, "TauIDSF_muToTau_up", w_m_up)
        self._set_tree_attr(tree, "TauIDSF_muToTau_down", w_m_dn)

        self._set_tree_attr(tree, "TauIDSF_otherVSlepWP_up", w_o_up)
        self._set_tree_attr(tree, "TauIDSF_otherVSlepWP_down", w_o_dn)

        self.out.fillBranch("TauIDSF", nominal)
        self.out.fillBranch("TauIDSF_Lepton1", float(sf1))
        self.out.fillBranch("TauIDSF_Lepton2", float(sf2))
        self.out.fillBranch("TauIDSF_Lepton3", float(sf3))

        self.out.fillBranch("TauIDSF_up", w_tot_up)
        self.out.fillBranch("TauIDSF_down", w_tot_dn)

        self.out.fillBranch("TauIDSF_genuineTau_up", w_g_up)
        self.out.fillBranch("TauIDSF_genuineTau_down", w_g_dn)

        self.out.fillBranch("TauIDSF_eToTau_up", w_e_up)
        self.out.fillBranch("TauIDSF_eToTau_down", w_e_dn)

        self.out.fillBranch("TauIDSF_muToTau_up", w_m_up)
        self.out.fillBranch("TauIDSF_muToTau_down", w_m_dn)

        self.out.fillBranch("TauIDSF_otherVSlepWP_up", w_o_up)
        self.out.fillBranch("TauIDSF_otherVSlepWP_down", w_o_dn)

        return True


def _factory(mode, year="2018"):
    return lambda: TauSFProducer(
        mode=mode,
        year=year,
        vsjet_wp="Medium",
        vsjet_wp_vse="VVLoose",
        vsjet_flag="pt",
        vse_wp="VVLoose",
        vsmu_wp="VLoose",
        apply_other_vslepwp_unc=True,
    )


MODES = [
    "ElElEl", "MuElEl", "ElMuMu", "MuMuMu",
    "TaTaTa", "ElTaTa", "MuTaTa", "TaElEl", "TaMuMu",
    "NPLElElEl", "NPLMuElEl", "NPLElMuMu", "NPLMuMuMu",
    "NPLTaTaTa", "NPLElTaTa", "NPLMuTaTa", "NPLTaElEl", "NPLTaMuMu",
]

for _mode in MODES:
    globals()["tauSF_%s_2018" % _mode] = _factory(_mode, year="2018")
