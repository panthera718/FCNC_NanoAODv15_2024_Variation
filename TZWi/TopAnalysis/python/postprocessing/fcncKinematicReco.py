import ROOT
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True

import os
from ROOT import TLorentzVector
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class FCNCKinematicReco(Module, object):
    def __init__(self, *args, **kwargs):
#        self.mode = kwargs.get("mode")
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch("KinTop_status", "i")

        self.out.branch("KinTopWb_pt", "F")
        self.out.branch("KinTopWb_eta", "F")
        self.out.branch("KinTopWb_phi", "F")
        self.out.branch("KinTopWb_mass", "F")

        self.out.branch("KinTopZq_pt", "F")
        self.out.branch("KinTopZq_eta", "F")
        self.out.branch("KinTopZq_phi", "F")
        self.out.branch("KinTopZq_mass", "F")
        self.out.branch("KinTop_deta", "F")
        self.out.branch("KinTop_dphi", "F")
        self.out.branch("KinTop_dR", "F")

        pass
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

    def analyze(self, event):

        # Constants
        m_W = 80.4  # W boson mass in GeV
        m_top = 172.5  # Top quark mass in GeV

        if event._tree.b_out_GoodLeptonCode == 111:

            Wlepvar = TLorentzVector()
            Zlep1var = TLorentzVector()
            Zlep2var = TLorentzVector()
            Wlepvar.SetPtEtaPhiM(event._tree.b_out_Lepton1_pt, event._tree.b_out_Lepton1_eta, event._tree.b_out_Lepton1_phi, event._tree.b_out_Lepton1_mass)
            Zlep1var.SetPtEtaPhiM(event._tree.b_out_Lepton2_pt, event._tree.b_out_Lepton2_eta, event._tree.b_out_Lepton2_phi, event._tree.b_out_Lepton2_mass)
            Zlep2var.SetPtEtaPhiM(event._tree.b_out_Lepton3_pt, event._tree.b_out_Lepton3_eta, event._tree.b_out_Lepton3_phi, event._tree.b_out_Lepton3_mass)

            setattr(event._tree, "b_out_Wlepvar", Wlepvar)
            setattr(event._tree, "b_out_Zlep1var", Zlep1var)
            setattr(event._tree, "b_out_Zlep2var", Zlep2var)

        if event._tree.b_out_GoodLeptonCode == 111 and event._tree.b_out_nGoodJet == 1 and event._tree.b_out_nBjet >= 1:

            jet_W = TLorentzVector()
            jet_W.SetPtEtaPhiM(event._tree.b_out_GoodJet_pt[0], event._tree.b_out_GoodJet_eta[0], event._tree.b_out_GoodJet_phi[0], event._tree.b_out_GoodJet_mass[0])
            jet_W_DeepFlavB = event._tree.b_out_GoodJet_DeepFlavB[0]
            jet_W_DeepFlavCvL = event._tree.b_out_GoodJet_DeepFlavCvL[0]
            jet_W_DeepFlavCvB = event._tree.b_out_GoodJet_DeepFlavCvB[0]
            jet_W_DeepFlavQG = event._tree.b_out_GoodJet_DeepFlavQG[0]

            wl_TLV = Wlepvar

            MET_pt = event._tree.b_out_MET_pt
            MET_phi = event._tree.b_out_MET_phi

            MET_px = MET_pt * math.cos(MET_phi)
            MET_py = MET_pt * math.sin(MET_phi)

            p_lepton = (Wlepvar.Px(), Wlepvar.Py(), Wlepvar.Pz(), Wlepvar.E())
            p_met = (MET_px, MET_py)
            p_nu_z_solutions = self.calculate_neutrino_pz(p_lepton, p_met, m_W)

            p_z_MET_plus = p_nu_z_solutions[0]
            p_z_MET_minus = p_nu_z_solutions[1]

            # Compute energies of MET
            E_MET_plus = math.sqrt(MET_pt ** 2 + p_z_MET_plus ** 2)
            E_MET_minus = math.sqrt(MET_pt ** 2 + p_z_MET_minus ** 2)

            # Create MET TLorentzVectors
            MET_TLV_plus = ROOT.TLorentzVector()
            MET_TLV_plus.SetPxPyPzE(MET_px, MET_py, p_z_MET_plus, E_MET_plus)

            MET_TLV_minus = ROOT.TLorentzVector()
            MET_TLV_minus.SetPxPyPzE(MET_px, MET_py, p_z_MET_minus, E_MET_minus)

            # Compute top_W_mass with MET_TLV_plus
            top_W_mass_plus = (jet_W + wl_TLV + MET_TLV_plus).M()
            # Compute top_W_mass with MET_TLV_minus
            top_W_mass_minus = (jet_W + wl_TLV + MET_TLV_minus).M()

            # Select the solution closer to m_top
            if abs(top_W_mass_plus - m_top) < abs(top_W_mass_minus - m_top):
                top_W_mass = top_W_mass_plus
                MET_TLV_used = MET_TLV_plus
            else:
                top_W_mass = top_W_mass_minus
                MET_TLV_used = MET_TLV_minus

            top_W_TLV = jet_W + wl_TLV + MET_TLV_used

            KinTopWb_pt   = top_W_TLV.Pt()
            KinTopWb_eta  = top_W_TLV.Eta()
            KinTopWb_phi  = top_W_TLV.Phi()
            KinTopWb_mass = top_W_TLV.M()

            setattr(event._tree, "b_out_jet_W", jet_W)
            setattr(event._tree, "b_out_jet_W_DeepFlavB", jet_W_DeepFlavB)
            setattr(event._tree, "b_out_jet_W_DeepFlavCvL", jet_W_DeepFlavCvL)
            setattr(event._tree, "b_out_jet_W_DeepFlavCvB", jet_W_DeepFlavCvB)
            setattr(event._tree, "b_out_jet_W_DeepFlavQG", jet_W_DeepFlavQG)
            setattr(event._tree, "b_out_jet_Z", 0)
            setattr(event._tree, "b_out_jet_Z_DeepFlavB", 0)
            setattr(event._tree, "b_out_jet_Z_DeepFlavCvL", 0)
            setattr(event._tree, "b_out_jet_Z_DeepFlavCvB", 0)
            setattr(event._tree, "b_out_jet_Z_DeepFlavQG", 0)

            setattr(event._tree, "b_out_KinTopWb_pt", KinTopWb_pt)
            setattr(event._tree, "b_out_KinTopWb_eta", KinTopWb_eta)
            setattr(event._tree, "b_out_KinTopWb_phi", KinTopWb_phi)
            setattr(event._tree, "b_out_KinTopWb_mass", KinTopWb_mass)

            setattr(event._tree, "b_out_KinTopZq_pt", 0)
            setattr(event._tree, "b_out_KinTopZq_eta", 0)
            setattr(event._tree, "b_out_KinTopZq_phi", 0)
            setattr(event._tree, "b_out_KinTopZq_mass", 0)
            setattr(event._tree, "b_out_KinTop_deta", 0)
            setattr(event._tree, "b_out_KinTop_dphi", 0)
            setattr(event._tree, "b_out_KinTop_dR", 0)

            self.out.fillBranch("KinTop_status", 1)

            self.out.fillBranch("KinTopWb_pt", KinTopWb_pt)
            self.out.fillBranch("KinTopWb_eta", KinTopWb_eta)
            self.out.fillBranch("KinTopWb_phi", KinTopWb_phi)
            self.out.fillBranch("KinTopWb_mass", KinTopWb_mass)

            self.out.fillBranch("KinTopZq_pt", 0)
            self.out.fillBranch("KinTopZq_eta", 0)
            self.out.fillBranch("KinTopZq_phi", 0)
            self.out.fillBranch("KinTopZq_mass", 0)
            self.out.fillBranch("KinTop_deta", 0)
            self.out.fillBranch("KinTop_dphi", 0)
            self.out.fillBranch("KinTop_dR", 0)

            return True

        elif event._tree.b_out_GoodLeptonCode == 111 and event._tree.b_out_nGoodJet >= 2 and event._tree.b_out_nBjet >= 1:

            # Given inputs (you need to provide actual TLorentzVectors and values)
            jet_TLV = []
            jet_DeepFlavB = []
            jet_DeepFlavCvL = []
            jet_DeepFlavCvB = []
            jet_DeepFlavQG = []
            for j_e in range(event._tree.b_out_nGoodJet):
                jet_TLV_element = TLorentzVector()
                jet_TLV_element.SetPtEtaPhiM(event._tree.b_out_GoodJet_pt[j_e], event._tree.b_out_GoodJet_eta[j_e], event._tree.b_out_GoodJet_phi[j_e], event._tree.b_out_GoodJet_mass[j_e])
                jet_TLV.append(jet_TLV_element)
                jet_DeepFlavB.append(event._tree.b_out_GoodJet_DeepFlavB[j_e])
                jet_DeepFlavCvL.append(event._tree.b_out_GoodJet_DeepFlavCvL[j_e])
                jet_DeepFlavCvB.append(event._tree.b_out_GoodJet_DeepFlavCvB[j_e])
                jet_DeepFlavQG.append(event._tree.b_out_GoodJet_DeepFlavQG[j_e])

            z_TLV = Zlep1var + Zlep2var
            wl_TLV = Wlepvar

            MET_pt = event._tree.b_out_MET_pt
            MET_phi = event._tree.b_out_MET_phi

            MET_px = MET_pt * math.cos(MET_phi)
            MET_py = MET_pt * math.sin(MET_phi)

            sigma_top_W = 1.
            sigma_top_Z = 1.

            p_lepton = (Wlepvar.Px(), Wlepvar.Py(), Wlepvar.Pz(), Wlepvar.E())
            p_met = (MET_px, MET_py)
            p_nu_z_solutions = self.calculate_neutrino_pz(p_lepton, p_met, m_W)

            p_z_MET_plus = p_nu_z_solutions[0]
            p_z_MET_minus = p_nu_z_solutions[1]

            # Compute energies of MET
            E_MET_plus = math.sqrt(MET_pt ** 2 + p_z_MET_plus ** 2)
            E_MET_minus = math.sqrt(MET_pt ** 2 + p_z_MET_minus ** 2)

            # Create MET TLorentzVectors
            MET_TLV_plus = ROOT.TLorentzVector()
            MET_TLV_plus.SetPxPyPzE(MET_px, MET_py, p_z_MET_plus, E_MET_plus)

            MET_TLV_minus = ROOT.TLorentzVector()
            MET_TLV_minus.SetPxPyPzE(MET_px, MET_py, p_z_MET_minus, E_MET_minus)

            # Initialize variables to store the minimum chi-square and corresponding indices
            min_chi_square = float('inf')
            best_i = -1
            best_j = -1
            best_MET_TLV = None  # To store the MET_TLV used for best_i and best_j
            best_top_W_mass = None
            best_top_Z_mass = None

            # Precompute top_Z_mass for each jet_TLV[j]
            top_Z_mass_list = []
            for j, jet_Z in enumerate(jet_TLV):
                top_Z_mass = (jet_Z + z_TLV).M()
                top_Z_mass_list.append(top_Z_mass)

            # Loop over jet_TLV[i] to compute top_W_mass
            for i, jet_W in enumerate(jet_TLV):

                # For each jet_TLV[j], compute chi-square
                for j, top_Z_mass in enumerate(top_Z_mass_list):

                    if i==j:
                        continue

                    # Compute top_W_mass with MET_TLV_plus
                    top_W_mass_plus = (jet_W + wl_TLV + MET_TLV_plus).M()
                    # Compute top_W_mass with MET_TLV_minus
                    top_W_mass_minus = (jet_W + wl_TLV + MET_TLV_minus).M()

                    # Select the solution closer to m_top
                    if abs(top_W_mass_plus - m_top) < abs(top_W_mass_minus - m_top):
                        top_W_mass = top_W_mass_plus
                        MET_TLV_used = MET_TLV_plus
                    else:
                        top_W_mass = top_W_mass_minus
                        MET_TLV_used = MET_TLV_minus

                    chi_square = ((m_top - top_W_mass) ** 2) / (sigma_top_W ** 2) + ((m_top - top_Z_mass) ** 2) / (sigma_top_Z ** 2)
                    # Update minimum chi-square and indices if current chi-square is lower
                    if chi_square < min_chi_square:
                        min_chi_square = chi_square
                        best_i = i
                        best_j = j
                        best_MET_TLV = MET_TLV_used
                        best_top_W_mass = top_W_mass
                        best_top_Z_mass = top_Z_mass

            top_W_TLV = jet_TLV[best_i] + wl_TLV + best_MET_TLV
            top_Z_TLV = jet_TLV[best_j] + z_TLV

            KinTopWb_pt   = top_W_TLV.Pt()
            KinTopWb_eta  = top_W_TLV.Eta()
            KinTopWb_phi  = top_W_TLV.Phi()
            KinTopWb_mass = top_W_TLV.M()
            KinTopZq_pt   = top_Z_TLV.Pt()
            KinTopZq_eta  = top_Z_TLV.Eta()
            KinTopZq_phi  = top_Z_TLV.Phi()
            KinTopZq_mass = top_Z_TLV.M()
            KinTop_deta   = abs(top_W_TLV.Eta()-top_Z_TLV.Eta()) 
            KinTop_dphi   = abs(top_W_TLV.Phi()-top_Z_TLV.Phi())
            KinTop_dR     = math.sqrt( (top_W_TLV.Eta()-top_Z_TLV.Eta())**2 + (top_W_TLV.Phi()-top_Z_TLV.Phi())**2 )

            setattr(event._tree, "b_out_jet_W", jet_TLV[best_i])
            setattr(event._tree, "b_out_jet_W_DeepFlavB", jet_DeepFlavB[best_i])
            setattr(event._tree, "b_out_jet_W_DeepFlavCvL", jet_DeepFlavCvL[best_i])
            setattr(event._tree, "b_out_jet_W_DeepFlavCvB", jet_DeepFlavCvB[best_i])
            setattr(event._tree, "b_out_jet_W_DeepFlavQG", jet_DeepFlavQG[best_i])
            setattr(event._tree, "b_out_jet_Z", jet_TLV[best_j])
            setattr(event._tree, "b_out_jet_Z_DeepFlavB", jet_DeepFlavB[best_j])
            setattr(event._tree, "b_out_jet_Z_DeepFlavCvL", jet_DeepFlavCvL[best_j])
            setattr(event._tree, "b_out_jet_Z_DeepFlavCvB", jet_DeepFlavCvB[best_j])
            setattr(event._tree, "b_out_jet_Z_DeepFlavQG", jet_DeepFlavQG[best_j])

            setattr(event._tree, "b_out_KinTopWb_pt", KinTopWb_pt)
            setattr(event._tree, "b_out_KinTopWb_eta", KinTopWb_eta)
            setattr(event._tree, "b_out_KinTopWb_phi", KinTopWb_phi)
            setattr(event._tree, "b_out_KinTopWb_mass", KinTopWb_mass)

            setattr(event._tree, "b_out_KinTopZq_pt", KinTopZq_pt)
            setattr(event._tree, "b_out_KinTopZq_eta", KinTopZq_eta)
            setattr(event._tree, "b_out_KinTopZq_phi", KinTopZq_phi)
            setattr(event._tree, "b_out_KinTopZq_mass", KinTopZq_mass)
            setattr(event._tree, "b_out_KinTop_deta", KinTop_deta)
            setattr(event._tree, "b_out_KinTop_dphi", KinTop_dphi)
            setattr(event._tree, "b_out_KinTop_dR", KinTop_dR)

            self.out.fillBranch("KinTop_status", 2)

            self.out.fillBranch("KinTopWb_pt", KinTopWb_pt)
            self.out.fillBranch("KinTopWb_eta", KinTopWb_eta)
            self.out.fillBranch("KinTopWb_phi", KinTopWb_phi)
            self.out.fillBranch("KinTopWb_mass", KinTopWb_mass)

            self.out.fillBranch("KinTopZq_pt", KinTopZq_pt)
            self.out.fillBranch("KinTopZq_eta", KinTopZq_eta)
            self.out.fillBranch("KinTopZq_phi", KinTopZq_phi)
            self.out.fillBranch("KinTopZq_mass", KinTopZq_mass)
            self.out.fillBranch("KinTop_deta", KinTop_deta)
            self.out.fillBranch("KinTop_dphi", KinTop_dphi)
            self.out.fillBranch("KinTop_dR", KinTop_dR)

            return True

        # Check basic event selection
        else:

            setattr(event._tree, "b_out_jet_W", 0)
            setattr(event._tree, "b_out_jet_W_DeepFlavB", 0)
            setattr(event._tree, "b_out_jet_Z", 0)
            setattr(event._tree, "b_out_jet_Z_DeepFlavB", 0)

            setattr(event._tree, "b_out_KinTopWb_pt", 0)
            setattr(event._tree, "b_out_KinTopWb_eta", 0)
            setattr(event._tree, "b_out_KinTopWb_phi", 0)
            setattr(event._tree, "b_out_KinTopWb_mass", 0)

            setattr(event._tree, "b_out_KinTopZq_pt", 0)
            setattr(event._tree, "b_out_KinTopZq_eta", 0)
            setattr(event._tree, "b_out_KinTopZq_phi", 0)
            setattr(event._tree, "b_out_KinTopZq_mass", 0)
            setattr(event._tree, "b_out_KinTop_deta", 0)
            setattr(event._tree, "b_out_KinTop_dphi", 0)
            setattr(event._tree, "b_out_KinTop_dR", 0)

            self.out.fillBranch("KinTop_status", 0)

            self.out.fillBranch("KinTopWb_pt", 0)
            self.out.fillBranch("KinTopWb_eta", 0)
            self.out.fillBranch("KinTopWb_phi", 0)
            self.out.fillBranch("KinTopWb_mass", 0)

            self.out.fillBranch("KinTopZq_pt", 0)
            self.out.fillBranch("KinTopZq_eta", 0)
            self.out.fillBranch("KinTopZq_phi", 0)
            self.out.fillBranch("KinTopZq_mass", 0)
            self.out.fillBranch("KinTop_deta", 0)
            self.out.fillBranch("KinTop_dphi", 0)
            self.out.fillBranch("KinTop_dR", 0)

            return True

fcncKinReco = lambda: FCNCKinematicReco()
