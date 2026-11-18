import ROOT
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True

import os
from ROOT import TLorentzVector
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class FCNCKinematicReco(Module, object):
    def __init__(self, *args, **kwargs):
        # Keep mode for interface compatibility with existing workflows.
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

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def _is_jetvar_mode(self):
        return self.jetVar is not None and self.jetSyst in ("Up", "Down")

    def _suffix(self):
        if not self._is_jetvar_mode():
            return ""
        return "_%s%s" % (self.jetVar, self.jetSyst)

    def _tree_name(self, base):
        return "b_out_%s%s" % (base, self._suffix())

    def _out_name(self, base):
        return "%s%s" % (base, self._suffix())

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch(self._out_name("KinTop_status"), "i")

        self.out.branch(self._out_name("KinTopWb_pt"), "F")
        self.out.branch(self._out_name("KinTopWb_eta"), "F")
        self.out.branch(self._out_name("KinTopWb_phi"), "F")
        self.out.branch(self._out_name("KinTopWb_mass"), "F")

        self.out.branch(self._out_name("KinTopZq_pt"), "F")
        self.out.branch(self._out_name("KinTopZq_eta"), "F")
        self.out.branch(self._out_name("KinTopZq_phi"), "F")
        self.out.branch(self._out_name("KinTopZq_mass"), "F")
        self.out.branch(self._out_name("KinTop_deta"), "F")
        self.out.branch(self._out_name("KinTop_dphi"), "F")
        self.out.branch(self._out_name("KinTop_dR"), "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        # Clean up attached TLorentzVectors on the PyROOT TTree proxy to prevent GC crashes
        if inputTree and hasattr(inputTree, 'tree') and inputTree.tree():
            tree = inputTree.tree()
            for attr in list(tree.__dict__.keys()):
                if attr.startswith("b_out_Wlepvar") or attr.startswith("b_out_Zlep") or attr.startswith("b_out_jet_"):
                    delattr(tree, attr)

    def calculate_neutrino_pz(self, p_lepton, p_met, m_W):
        lep_px, lep_py, lep_pz, lep_E = p_lepton
        met_px, met_py = p_met

        lep_pt_sq = lep_px**2 + lep_py**2
        met_pt_sq = met_px**2 + met_py**2
        lep_dot_met = lep_px*met_px + lep_py*met_py
        det_in = (m_W*m_W + 2*lep_dot_met)**2 - 4*lep_pt_sq*met_pt_sq
        det = max(0., det_in)

        nu_pz1 = ((m_W*m_W + 2*lep_dot_met)*lep_pz + lep_E*math.sqrt(det))/(2*lep_pt_sq)
        nu_pz2 = ((m_W*m_W + 2*lep_dot_met)*lep_pz - lep_E*math.sqrt(det))/(2*lep_pt_sq)

        return [nu_pz1, nu_pz2]

    def _get_value(self, tree, base):
        return getattr(tree, self._tree_name(base))

    def _set_value(self, tree, base, value):
        setattr(tree, self._tree_name(base), value)

    def _fill_output(self, status, KinTopWb_pt, KinTopWb_eta, KinTopWb_phi, KinTopWb_mass,
                     KinTopZq_pt, KinTopZq_eta, KinTopZq_phi, KinTopZq_mass,
                     KinTop_deta, KinTop_dphi, KinTop_dR):
        self.out.fillBranch(self._out_name("KinTop_status"), status)

        self.out.fillBranch(self._out_name("KinTopWb_pt"), KinTopWb_pt)
        self.out.fillBranch(self._out_name("KinTopWb_eta"), KinTopWb_eta)
        self.out.fillBranch(self._out_name("KinTopWb_phi"), KinTopWb_phi)
        self.out.fillBranch(self._out_name("KinTopWb_mass"), KinTopWb_mass)

        self.out.fillBranch(self._out_name("KinTopZq_pt"), KinTopZq_pt)
        self.out.fillBranch(self._out_name("KinTopZq_eta"), KinTopZq_eta)
        self.out.fillBranch(self._out_name("KinTopZq_phi"), KinTopZq_phi)
        self.out.fillBranch(self._out_name("KinTopZq_mass"), KinTopZq_mass)
        self.out.fillBranch(self._out_name("KinTop_deta"), KinTop_deta)
        self.out.fillBranch(self._out_name("KinTop_dphi"), KinTop_dphi)
        self.out.fillBranch(self._out_name("KinTop_dR"), KinTop_dR)

    def analyze(self, event):
        # Constants
        m_W = 80.4
        m_top = 172.5
        tree = event._tree

        if tree.b_out_GoodLeptonCode == 111:
            Wlepvar = TLorentzVector()
            Zlep1var = TLorentzVector()
            Zlep2var = TLorentzVector()
            Wlepvar.SetPtEtaPhiM(tree.b_out_Lepton1_pt, tree.b_out_Lepton1_eta, tree.b_out_Lepton1_phi, tree.b_out_Lepton1_mass)
            Zlep1var.SetPtEtaPhiM(tree.b_out_Lepton2_pt, tree.b_out_Lepton2_eta, tree.b_out_Lepton2_phi, tree.b_out_Lepton2_mass)
            Zlep2var.SetPtEtaPhiM(tree.b_out_Lepton3_pt, tree.b_out_Lepton3_eta, tree.b_out_Lepton3_phi, tree.b_out_Lepton3_mass)

            self._set_value(tree, "Wlepvar", Wlepvar)
            self._set_value(tree, "Zlep1var", Zlep1var)
            self._set_value(tree, "Zlep2var", Zlep2var)

        nGoodJet = self._get_value(tree, "nGoodJet") if self._is_jetvar_mode() else tree.b_out_nGoodJet
        nBjet = self._get_value(tree, "nBjet") if self._is_jetvar_mode() else tree.b_out_nBjet

        if tree.b_out_GoodLeptonCode == 111 and nGoodJet == 1 and nBjet >= 1:
            jet_W = TLorentzVector()
            jet_W.SetPtEtaPhiM(
                self._get_value(tree, "GoodJet_pt")[0],
                self._get_value(tree, "GoodJet_eta")[0],
                self._get_value(tree, "GoodJet_phi")[0],
                self._get_value(tree, "GoodJet_mass")[0],
            ) if self._is_jetvar_mode() else jet_W.SetPtEtaPhiM(
                tree.b_out_GoodJet_pt[0],
                tree.b_out_GoodJet_eta[0],
                tree.b_out_GoodJet_phi[0],
                tree.b_out_GoodJet_mass[0],
            )

            if self._is_jetvar_mode():
                jet_W_DeepFlavB = self._get_value(tree, "GoodJet_DeepFlavB")[0]
                jet_W_DeepFlavCvL = self._get_value(tree, "GoodJet_DeepFlavCvL")[0]
                jet_W_DeepFlavCvB = self._get_value(tree, "GoodJet_DeepFlavCvB")[0]
                jet_W_DeepFlavQG = self._get_value(tree, "GoodJet_DeepFlavQG")[0]
                MET_pt = self._get_value(tree, "MET_pt")
                MET_phi = self._get_value(tree, "MET_phi")
            else:
                jet_W_DeepFlavB = tree.b_out_GoodJet_DeepFlavB[0]
                jet_W_DeepFlavCvL = tree.b_out_GoodJet_DeepFlavCvL[0]
                jet_W_DeepFlavCvB = tree.b_out_GoodJet_DeepFlavCvB[0]
                jet_W_DeepFlavQG = tree.b_out_GoodJet_DeepFlavQG[0]
                MET_pt = tree.b_out_MET_pt
                MET_phi = tree.b_out_MET_phi

            wl_TLV = Wlepvar

            MET_px = MET_pt * math.cos(MET_phi)
            MET_py = MET_pt * math.sin(MET_phi)

            p_lepton = (Wlepvar.Px(), Wlepvar.Py(), Wlepvar.Pz(), Wlepvar.E())
            p_met = (MET_px, MET_py)
            p_nu_z_solutions = self.calculate_neutrino_pz(p_lepton, p_met, m_W)

            p_z_MET_plus = p_nu_z_solutions[0]
            p_z_MET_minus = p_nu_z_solutions[1]

            E_MET_plus = math.sqrt(MET_pt ** 2 + p_z_MET_plus ** 2)
            E_MET_minus = math.sqrt(MET_pt ** 2 + p_z_MET_minus ** 2)

            MET_TLV_plus = ROOT.TLorentzVector()
            MET_TLV_plus.SetPxPyPzE(MET_px, MET_py, p_z_MET_plus, E_MET_plus)

            MET_TLV_minus = ROOT.TLorentzVector()
            MET_TLV_minus.SetPxPyPzE(MET_px, MET_py, p_z_MET_minus, E_MET_minus)

            top_W_mass_plus = (jet_W + wl_TLV + MET_TLV_plus).M()
            top_W_mass_minus = (jet_W + wl_TLV + MET_TLV_minus).M()

            if abs(top_W_mass_plus - m_top) < abs(top_W_mass_minus - m_top):
                MET_TLV_used = MET_TLV_plus
            else:
                MET_TLV_used = MET_TLV_minus

            top_W_TLV = jet_W + wl_TLV + MET_TLV_used

            KinTopWb_pt = top_W_TLV.Pt()
            KinTopWb_eta = top_W_TLV.Eta()
            KinTopWb_phi = top_W_TLV.Phi()
            KinTopWb_mass = top_W_TLV.M()

            self._set_value(tree, "jet_W", jet_W)
            self._set_value(tree, "jet_W_DeepFlavB", jet_W_DeepFlavB)
            self._set_value(tree, "jet_W_DeepFlavCvL", jet_W_DeepFlavCvL)
            self._set_value(tree, "jet_W_DeepFlavCvB", jet_W_DeepFlavCvB)
            self._set_value(tree, "jet_W_DeepFlavQG", jet_W_DeepFlavQG)
            self._set_value(tree, "jet_Z", 0)
            self._set_value(tree, "jet_Z_DeepFlavB", 0)
            self._set_value(tree, "jet_Z_DeepFlavCvL", 0)
            self._set_value(tree, "jet_Z_DeepFlavCvB", 0)
            self._set_value(tree, "jet_Z_DeepFlavQG", 0)

            self._set_value(tree, "KinTopWb_pt", KinTopWb_pt)
            self._set_value(tree, "KinTopWb_eta", KinTopWb_eta)
            self._set_value(tree, "KinTopWb_phi", KinTopWb_phi)
            self._set_value(tree, "KinTopWb_mass", KinTopWb_mass)

            self._set_value(tree, "KinTopZq_pt", 0)
            self._set_value(tree, "KinTopZq_eta", 0)
            self._set_value(tree, "KinTopZq_phi", 0)
            self._set_value(tree, "KinTopZq_mass", 0)
            self._set_value(tree, "KinTop_deta", 0)
            self._set_value(tree, "KinTop_dphi", 0)
            self._set_value(tree, "KinTop_dR", 0)

            self._fill_output(1, KinTopWb_pt, KinTopWb_eta, KinTopWb_phi, KinTopWb_mass, 0, 0, 0, 0, 0, 0, 0)
            return True

        elif tree.b_out_GoodLeptonCode == 111 and nGoodJet >= 2 and nBjet >= 1:
            jet_TLV = []
            jet_DeepFlavB = []
            jet_DeepFlavCvL = []
            jet_DeepFlavCvB = []
            jet_DeepFlavQG = []

            if self._is_jetvar_mode():
                goodjet_pt = self._get_value(tree, "GoodJet_pt")
                goodjet_eta = self._get_value(tree, "GoodJet_eta")
                goodjet_phi = self._get_value(tree, "GoodJet_phi")
                goodjet_mass = self._get_value(tree, "GoodJet_mass")
                goodjet_b = self._get_value(tree, "GoodJet_DeepFlavB")
                goodjet_cvl = self._get_value(tree, "GoodJet_DeepFlavCvL")
                goodjet_cvb = self._get_value(tree, "GoodJet_DeepFlavCvB")
                goodjet_qg = self._get_value(tree, "GoodJet_DeepFlavQG")
                MET_pt = self._get_value(tree, "MET_pt")
                MET_phi = self._get_value(tree, "MET_phi")
            else:
                goodjet_pt = tree.b_out_GoodJet_pt
                goodjet_eta = tree.b_out_GoodJet_eta
                goodjet_phi = tree.b_out_GoodJet_phi
                goodjet_mass = tree.b_out_GoodJet_mass
                goodjet_b = tree.b_out_GoodJet_DeepFlavB
                goodjet_cvl = tree.b_out_GoodJet_DeepFlavCvL
                goodjet_cvb = tree.b_out_GoodJet_DeepFlavCvB
                goodjet_qg = tree.b_out_GoodJet_DeepFlavQG
                MET_pt = tree.b_out_MET_pt
                MET_phi = tree.b_out_MET_phi

            for j_e in range(nGoodJet):
                jet_TLV_element = TLorentzVector()
                jet_TLV_element.SetPtEtaPhiM(goodjet_pt[j_e], goodjet_eta[j_e], goodjet_phi[j_e], goodjet_mass[j_e])
                jet_TLV.append(jet_TLV_element)
                jet_DeepFlavB.append(goodjet_b[j_e])
                jet_DeepFlavCvL.append(goodjet_cvl[j_e])
                jet_DeepFlavCvB.append(goodjet_cvb[j_e])
                jet_DeepFlavQG.append(goodjet_qg[j_e])

            z_TLV = Zlep1var + Zlep2var
            wl_TLV = Wlepvar

            MET_px = MET_pt * math.cos(MET_phi)
            MET_py = MET_pt * math.sin(MET_phi)

            sigma_top_W = 1.
            sigma_top_Z = 1.

            p_lepton = (Wlepvar.Px(), Wlepvar.Py(), Wlepvar.Pz(), Wlepvar.E())
            p_met = (MET_px, MET_py)
            p_nu_z_solutions = self.calculate_neutrino_pz(p_lepton, p_met, m_W)

            p_z_MET_plus = p_nu_z_solutions[0]
            p_z_MET_minus = p_nu_z_solutions[1]

            E_MET_plus = math.sqrt(MET_pt ** 2 + p_z_MET_plus ** 2)
            E_MET_minus = math.sqrt(MET_pt ** 2 + p_z_MET_minus ** 2)

            MET_TLV_plus = ROOT.TLorentzVector()
            MET_TLV_plus.SetPxPyPzE(MET_px, MET_py, p_z_MET_plus, E_MET_plus)

            MET_TLV_minus = ROOT.TLorentzVector()
            MET_TLV_minus.SetPxPyPzE(MET_px, MET_py, p_z_MET_minus, E_MET_minus)

            min_chi_square = float('inf')
            best_i = -1
            best_j = -1
            best_MET_TLV = None

            top_Z_mass_list = []
            for j, jet_Z in enumerate(jet_TLV):
                top_Z_mass = (jet_Z + z_TLV).M()
                top_Z_mass_list.append(top_Z_mass)

            for i, jet_W in enumerate(jet_TLV):
                for j, top_Z_mass in enumerate(top_Z_mass_list):
                    if i == j:
                        continue

                    top_W_mass_plus = (jet_W + wl_TLV + MET_TLV_plus).M()
                    top_W_mass_minus = (jet_W + wl_TLV + MET_TLV_minus).M()

                    if abs(top_W_mass_plus - m_top) < abs(top_W_mass_minus - m_top):
                        top_W_mass = top_W_mass_plus
                        MET_TLV_used = MET_TLV_plus
                    else:
                        top_W_mass = top_W_mass_minus
                        MET_TLV_used = MET_TLV_minus

                    chi_square = ((m_top - top_W_mass) ** 2) / (sigma_top_W ** 2) + ((m_top - top_Z_mass) ** 2) / (sigma_top_Z ** 2)
                    if chi_square < min_chi_square:
                        min_chi_square = chi_square
                        best_i = i
                        best_j = j
                        best_MET_TLV = MET_TLV_used

            top_W_TLV = jet_TLV[best_i] + wl_TLV + best_MET_TLV
            top_Z_TLV = jet_TLV[best_j] + z_TLV

            KinTopWb_pt = top_W_TLV.Pt()
            KinTopWb_eta = top_W_TLV.Eta()
            KinTopWb_phi = top_W_TLV.Phi()
            KinTopWb_mass = top_W_TLV.M()
            KinTopZq_pt = top_Z_TLV.Pt()
            KinTopZq_eta = top_Z_TLV.Eta()
            KinTopZq_phi = top_Z_TLV.Phi()
            KinTopZq_mass = top_Z_TLV.M()
            KinTop_deta = abs(top_W_TLV.Eta() - top_Z_TLV.Eta())
            KinTop_dphi = abs(top_W_TLV.Phi() - top_Z_TLV.Phi())
            KinTop_dR = math.sqrt((top_W_TLV.Eta() - top_Z_TLV.Eta())**2 + (top_W_TLV.Phi() - top_Z_TLV.Phi())**2)

            self._set_value(tree, "jet_W", jet_TLV[best_i])
            self._set_value(tree, "jet_W_DeepFlavB", jet_DeepFlavB[best_i])
            self._set_value(tree, "jet_W_DeepFlavCvL", jet_DeepFlavCvL[best_i])
            self._set_value(tree, "jet_W_DeepFlavCvB", jet_DeepFlavCvB[best_i])
            self._set_value(tree, "jet_W_DeepFlavQG", jet_DeepFlavQG[best_i])
            self._set_value(tree, "jet_Z", jet_TLV[best_j])
            self._set_value(tree, "jet_Z_DeepFlavB", jet_DeepFlavB[best_j])
            self._set_value(tree, "jet_Z_DeepFlavCvL", jet_DeepFlavCvL[best_j])
            self._set_value(tree, "jet_Z_DeepFlavCvB", jet_DeepFlavCvB[best_j])
            self._set_value(tree, "jet_Z_DeepFlavQG", jet_DeepFlavQG[best_j])

            self._set_value(tree, "KinTopWb_pt", KinTopWb_pt)
            self._set_value(tree, "KinTopWb_eta", KinTopWb_eta)
            self._set_value(tree, "KinTopWb_phi", KinTopWb_phi)
            self._set_value(tree, "KinTopWb_mass", KinTopWb_mass)

            self._set_value(tree, "KinTopZq_pt", KinTopZq_pt)
            self._set_value(tree, "KinTopZq_eta", KinTopZq_eta)
            self._set_value(tree, "KinTopZq_phi", KinTopZq_phi)
            self._set_value(tree, "KinTopZq_mass", KinTopZq_mass)
            self._set_value(tree, "KinTop_deta", KinTop_deta)
            self._set_value(tree, "KinTop_dphi", KinTop_dphi)
            self._set_value(tree, "KinTop_dR", KinTop_dR)

            self._fill_output(2, KinTopWb_pt, KinTopWb_eta, KinTopWb_phi, KinTopWb_mass,
                              KinTopZq_pt, KinTopZq_eta, KinTopZq_phi, KinTopZq_mass,
                              KinTop_deta, KinTop_dphi, KinTop_dR)
            return True

        else:
            self._set_value(tree, "jet_W", 0)
            self._set_value(tree, "jet_W_DeepFlavB", 0)
            self._set_value(tree, "jet_Z", 0)
            self._set_value(tree, "jet_Z_DeepFlavB", 0)

            self._set_value(tree, "KinTopWb_pt", 0)
            self._set_value(tree, "KinTopWb_eta", 0)
            self._set_value(tree, "KinTopWb_phi", 0)
            self._set_value(tree, "KinTopWb_mass", 0)

            self._set_value(tree, "KinTopZq_pt", 0)
            self._set_value(tree, "KinTopZq_eta", 0)
            self._set_value(tree, "KinTopZq_phi", 0)
            self._set_value(tree, "KinTopZq_mass", 0)
            self._set_value(tree, "KinTop_deta", 0)
            self._set_value(tree, "KinTop_dphi", 0)
            self._set_value(tree, "KinTop_dR", 0)

            self._fill_output(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
            return True


fcncKinReco = lambda: FCNCKinematicReco()

JETVARS = [
    "AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias",
    "Fragmentation", "SinglePionECAL", "SinglePionHCAL",
    "FlavorQCD", "TimePtEta",
    "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
    "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF",
    "RelativeBal", "RelativeSample", "RelativeFSR",
    "RelativeStatFSR", "RelativeStatEC", "RelativeStatHF",
    "PileUpDataMC", "PileUpPtRef", "PileUpPtBB",
    "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF",
]

for _jetvar in JETVARS:
    globals()["fcncKinRecoUp_%s" % _jetvar] = (lambda jv=_jetvar: FCNCKinematicReco(jetVar=jv, jetSyst="Up"))
    globals()["fcncKinRecoDown_%s" % _jetvar] = (lambda jv=_jetvar: FCNCKinematicReco(jetVar=jv, jetSyst="Down"))
