import ROOT
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True

import os
from ROOT import TLorentzVector
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class FCNCKinematicReco(Module, object):
    def __init__(self, *args, **kwargs):
        self.mode = kwargs.get("mode")
        self.jetVar = kwargs.get("jetVar")
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        jetVar = self.jetVar

        self.out.branch("KinTop_status_%sUp" % jetVar, "i")

        self.out.branch("KinTopWb_pt_%sUp" % jetVar, "F")
        self.out.branch("KinTopWb_eta_%sUp" % jetVar, "F")
        self.out.branch("KinTopWb_phi_%sUp" % jetVar, "F")
        self.out.branch("KinTopWb_mass_%sUp" % jetVar, "F")

        self.out.branch("KinTopZq_pt_%sUp" % jetVar, "F")
        self.out.branch("KinTopZq_eta_%sUp" % jetVar, "F")
        self.out.branch("KinTopZq_phi_%sUp" % jetVar, "F")
        self.out.branch("KinTopZq_mass_%sUp" % jetVar, "F")
        self.out.branch("KinTop_deta_%sUp" % jetVar, "F")
        self.out.branch("KinTop_dphi_%sUp" % jetVar, "F")
        self.out.branch("KinTop_dR_%sUp" % jetVar, "F")

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
        jetVar = self.jetVar

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

            setattr(event._tree, "b_out_Wlepvar_%sUp" % jetVar, Wlepvar)
            setattr(event._tree, "b_out_Zlep1var_%sUp" % jetVar, Zlep1var)
            setattr(event._tree, "b_out_Zlep2var_%sUp" % jetVar, Zlep2var)

        if event._tree.b_out_GoodLeptonCode == 111 and getattr(event._tree, "b_out_nGoodJet_%sUp" % jetVar) == 1 and getattr(event._tree, "b_out_nBjet_%sUp" % jetVar) >= 1:

            jet_W = TLorentzVector()
            jet_W.SetPtEtaPhiM(getattr(event._tree, "b_out_GoodJet_pt_%sUp" % jetVar)[0], getattr(event._tree, "b_out_GoodJet_eta_%sUp" % jetVar)[0], getattr(event._tree, "b_out_GoodJet_phi_%sUp" % jetVar)[0], getattr(event._tree, "b_out_GoodJet_mass_%sUp" % jetVar)[0])
            jet_W_DeepFlavB = getattr(event._tree, "b_out_GoodJet_DeepFlavB_%sUp" % jetVar)[0]
            jet_W_DeepFlavCvL = getattr(event._tree, "b_out_GoodJet_DeepFlavCvL_%sUp" % jetVar)[0]
            jet_W_DeepFlavCvB = getattr(event._tree, "b_out_GoodJet_DeepFlavCvB_%sUp" % jetVar)[0]
            jet_W_DeepFlavQG = getattr(event._tree, "b_out_GoodJet_DeepFlavQG_%sUp" % jetVar)[0]

            wl_TLV = Wlepvar

            MET_pt = getattr(event._tree, "b_out_MET_pt_%sUp" % jetVar)
            MET_phi = getattr(event._tree, "b_out_MET_phi_%sUp" % jetVar)

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

            setattr(event._tree, "b_out_jet_W_%sUp" % jetVar, jet_W)
            setattr(event._tree, "b_out_jet_W_DeepFlavB_%sUp" % jetVar, jet_W_DeepFlavB)
            setattr(event._tree, "b_out_jet_W_DeepFlavCvL_%sUp" % jetVar, jet_W_DeepFlavCvL)
            setattr(event._tree, "b_out_jet_W_DeepFlavCvB_%sUp" % jetVar, jet_W_DeepFlavCvB)
            setattr(event._tree, "b_out_jet_W_DeepFlavQG_%sUp" % jetVar, jet_W_DeepFlavQG)
            setattr(event._tree, "b_out_jet_Z_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_jet_Z_DeepFlavB_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_jet_Z_DeepFlavCvL_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_jet_Z_DeepFlavCvB_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_jet_Z_DeepFlavQG_%sUp" % jetVar, 0)

            setattr(event._tree, "b_out_KinTopWb_pt_%sUp" % jetVar, KinTopWb_pt)
            setattr(event._tree, "b_out_KinTopWb_eta_%sUp" % jetVar, KinTopWb_eta)
            setattr(event._tree, "b_out_KinTopWb_phi_%sUp" % jetVar, KinTopWb_phi)
            setattr(event._tree, "b_out_KinTopWb_mass_%sUp" % jetVar, KinTopWb_mass)

            setattr(event._tree, "b_out_KinTopZq_pt_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_KinTopZq_eta_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_KinTopZq_phi_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_KinTopZq_mass_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_KinTop_deta_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_KinTop_dphi_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_KinTop_dR_%sUp" % jetVar, 0)

            self.out.fillBranch("KinTop_status_%sUp" % jetVar, 1)

            self.out.fillBranch("KinTopWb_pt_%sUp" % jetVar, KinTopWb_pt)
            self.out.fillBranch("KinTopWb_eta_%sUp" % jetVar, KinTopWb_eta)
            self.out.fillBranch("KinTopWb_phi_%sUp" % jetVar, KinTopWb_phi)
            self.out.fillBranch("KinTopWb_mass_%sUp" % jetVar, KinTopWb_mass)

            self.out.fillBranch("KinTopZq_pt_%sUp" % jetVar, 0)
            self.out.fillBranch("KinTopZq_eta_%sUp" % jetVar, 0)
            self.out.fillBranch("KinTopZq_phi_%sUp" % jetVar, 0)
            self.out.fillBranch("KinTopZq_mass_%sUp" % jetVar, 0)
            self.out.fillBranch("KinTop_deta_%sUp" % jetVar, 0)
            self.out.fillBranch("KinTop_dphi_%sUp" % jetVar, 0)
            self.out.fillBranch("KinTop_dR_%sUp" % jetVar, 0)

            return True

        elif event._tree.b_out_GoodLeptonCode == 111 and getattr(event._tree, "b_out_nGoodJet_%sUp" % jetVar) >= 2 and getattr(event._tree, "b_out_nBjet_%sUp" % jetVar) >= 1:

            jet_TLV = []
            jet_DeepFlavB = []
            jet_DeepFlavCvL = []
            jet_DeepFlavCvB = []
            jet_DeepFlavQG = []
            for j_e in range(getattr(event._tree, "b_out_nGoodJet_%sUp" % jetVar)):
                jet_TLV_element = TLorentzVector()
                jet_TLV_element.SetPtEtaPhiM(getattr(event._tree, "b_out_GoodJet_pt_%sUp" % jetVar)[j_e], getattr(event._tree, "b_out_GoodJet_eta_%sUp" % jetVar)[j_e], getattr(event._tree, "b_out_GoodJet_phi_%sUp" % jetVar)[j_e], getattr(event._tree, "b_out_GoodJet_mass_%sUp" % jetVar)[j_e])
                jet_TLV.append(jet_TLV_element)
                jet_DeepFlavB.append(getattr(event._tree, "b_out_GoodJet_DeepFlavB_%sUp" % jetVar)[j_e])
                jet_DeepFlavCvL.append(getattr(event._tree, "b_out_GoodJet_DeepFlavCvL_%sUp" % jetVar)[j_e])
                jet_DeepFlavCvB.append(getattr(event._tree, "b_out_GoodJet_DeepFlavCvB_%sUp" % jetVar)[j_e])
                jet_DeepFlavQG.append(getattr(event._tree, "b_out_GoodJet_DeepFlavQG_%sUp" % jetVar)[j_e])

            z_TLV = Zlep1var + Zlep2var
            wl_TLV = Wlepvar

            MET_pt = getattr(event._tree, "b_out_MET_pt_%sUp" % jetVar)
            MET_phi = getattr(event._tree, "b_out_MET_phi_%sUp" % jetVar)

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

            setattr(event._tree, "b_out_jet_W_%sUp" % jetVar, jet_TLV[best_i])
            setattr(event._tree, "b_out_jet_W_DeepFlavB_%sUp" % jetVar, jet_DeepFlavB[best_i])
            setattr(event._tree, "b_out_jet_W_DeepFlavCvL_%sUp" % jetVar, jet_DeepFlavCvL[best_i])
            setattr(event._tree, "b_out_jet_W_DeepFlavCvB_%sUp" % jetVar, jet_DeepFlavCvB[best_i])
            setattr(event._tree, "b_out_jet_W_DeepFlavQG_%sUp" % jetVar, jet_DeepFlavQG[best_i])
            setattr(event._tree, "b_out_jet_Z_%sUp" % jetVar, jet_TLV[best_j])
            setattr(event._tree, "b_out_jet_Z_DeepFlavB_%sUp" % jetVar, jet_DeepFlavB[best_j])
            setattr(event._tree, "b_out_jet_Z_DeepFlavCvL_%sUp" % jetVar, jet_DeepFlavCvL[best_j])
            setattr(event._tree, "b_out_jet_Z_DeepFlavCvB_%sUp" % jetVar, jet_DeepFlavCvB[best_j])
            setattr(event._tree, "b_out_jet_Z_DeepFlavQG_%sUp" % jetVar, jet_DeepFlavQG[best_j])

            setattr(event._tree, "b_out_KinTopWb_pt_%sUp" % jetVar, KinTopWb_pt)
            setattr(event._tree, "b_out_KinTopWb_eta_%sUp" % jetVar, KinTopWb_eta)
            setattr(event._tree, "b_out_KinTopWb_phi_%sUp" % jetVar, KinTopWb_phi)
            setattr(event._tree, "b_out_KinTopWb_mass_%sUp" % jetVar, KinTopWb_mass)

            setattr(event._tree, "b_out_KinTopZq_pt_%sUp" % jetVar, KinTopZq_pt)
            setattr(event._tree, "b_out_KinTopZq_eta_%sUp" % jetVar, KinTopZq_eta)
            setattr(event._tree, "b_out_KinTopZq_phi_%sUp" % jetVar, KinTopZq_phi)
            setattr(event._tree, "b_out_KinTopZq_mass_%sUp" % jetVar, KinTopZq_mass)
            setattr(event._tree, "b_out_KinTop_deta_%sUp" % jetVar, KinTop_deta)
            setattr(event._tree, "b_out_KinTop_dphi_%sUp" % jetVar, KinTop_dphi)
            setattr(event._tree, "b_out_KinTop_dR_%sUp" % jetVar, KinTop_dR)

            self.out.fillBranch("KinTop_status_%sUp" % jetVar, 2)

            self.out.fillBranch("KinTopWb_pt_%sUp" % jetVar, KinTopWb_pt)
            self.out.fillBranch("KinTopWb_eta_%sUp" % jetVar, KinTopWb_eta)
            self.out.fillBranch("KinTopWb_phi_%sUp" % jetVar, KinTopWb_phi)
            self.out.fillBranch("KinTopWb_mass_%sUp" % jetVar, KinTopWb_mass)

            self.out.fillBranch("KinTopZq_pt_%sUp" % jetVar, KinTopZq_pt)
            self.out.fillBranch("KinTopZq_eta_%sUp" % jetVar, KinTopZq_eta)
            self.out.fillBranch("KinTopZq_phi_%sUp" % jetVar, KinTopZq_phi)
            self.out.fillBranch("KinTopZq_mass_%sUp" % jetVar, KinTopZq_mass)
            self.out.fillBranch("KinTop_deta_%sUp" % jetVar, KinTop_deta)
            self.out.fillBranch("KinTop_dphi_%sUp" % jetVar, KinTop_dphi)
            self.out.fillBranch("KinTop_dR_%sUp" % jetVar, KinTop_dR)

            return True

        # Check basic event selection
        else:

            setattr(event._tree, "b_out_jet_W_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_jet_W_DeepFlavB_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_jet_Z_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_jet_Z_DeepFlavB_%sUp" % jetVar, 0)

            setattr(event._tree, "b_out_KinTopWb_pt_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_KinTopWb_eta_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_KinTopWb_phi_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_KinTopWb_mass_%sUp" % jetVar, 0)

            setattr(event._tree, "b_out_KinTopZq_pt_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_KinTopZq_eta_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_KinTopZq_phi_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_KinTopZq_mass_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_KinTop_deta_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_KinTop_dphi_%sUp" % jetVar, 0)
            setattr(event._tree, "b_out_KinTop_dR_%sUp" % jetVar, 0)

            self.out.fillBranch("KinTop_status_%sUp" % jetVar, 0)

            self.out.fillBranch("KinTopWb_pt_%sUp" % jetVar, 0)
            self.out.fillBranch("KinTopWb_eta_%sUp" % jetVar, 0)
            self.out.fillBranch("KinTopWb_phi_%sUp" % jetVar, 0)
            self.out.fillBranch("KinTopWb_mass_%sUp" % jetVar, 0)

            self.out.fillBranch("KinTopZq_pt_%sUp" % jetVar, 0)
            self.out.fillBranch("KinTopZq_eta_%sUp" % jetVar, 0)
            self.out.fillBranch("KinTopZq_phi_%sUp" % jetVar, 0)
            self.out.fillBranch("KinTopZq_mass_%sUp" % jetVar, 0)
            self.out.fillBranch("KinTop_deta_%sUp" % jetVar, 0)
            self.out.fillBranch("KinTop_dphi_%sUp" % jetVar, 0)
            self.out.fillBranch("KinTop_dR_%sUp" % jetVar, 0)

            return True

fcncKinReco_AbsoluteStat = lambda : FCNCKinematicReco(jetVar="AbsoluteStat")
fcncKinReco_AbsoluteScale = lambda : FCNCKinematicReco(jetVar="AbsoluteScale")
fcncKinReco_AbsoluteMPFBias = lambda : FCNCKinematicReco(jetVar="AbsoluteMPFBias")
fcncKinReco_Fragmentation = lambda : FCNCKinematicReco(jetVar="Fragmentation")
fcncKinReco_SinglePionECAL = lambda : FCNCKinematicReco(jetVar="SinglePionECAL")
fcncKinReco_SinglePionHCAL = lambda : FCNCKinematicReco(jetVar="SinglePionHCAL")
fcncKinReco_FlavorQCD = lambda : FCNCKinematicReco(jetVar="FlavorQCD")
fcncKinReco_TimePtEta = lambda : FCNCKinematicReco(jetVar="TimePtEta")
fcncKinReco_RelativeJEREC1 = lambda : FCNCKinematicReco(jetVar="RelativeJEREC1")
fcncKinReco_RelativeJEREC2 = lambda : FCNCKinematicReco(jetVar="RelativeJEREC2")
fcncKinReco_RelativeJERHF = lambda : FCNCKinematicReco(jetVar="RelativeJERHF")
fcncKinReco_RelativePtBB = lambda : FCNCKinematicReco(jetVar="RelativePtBB")
fcncKinReco_RelativePtEC1 = lambda : FCNCKinematicReco(jetVar="RelativePtEC1")
fcncKinReco_RelativePtEC2 = lambda : FCNCKinematicReco(jetVar="RelativePtEC2")
fcncKinReco_RelativePtHF = lambda : FCNCKinematicReco(jetVar="RelativePtHF")
fcncKinReco_RelativeBal = lambda : FCNCKinematicReco(jetVar="RelativeBal")
fcncKinReco_RelativeSample = lambda : FCNCKinematicReco(jetVar="RelativeSample")
fcncKinReco_RelativeFSR = lambda : FCNCKinematicReco(jetVar="RelativeFSR")
fcncKinReco_RelativeStatFSR = lambda : FCNCKinematicReco(jetVar="RelativeStatFSR")
fcncKinReco_RelativeStatEC = lambda : FCNCKinematicReco(jetVar="RelativeStatEC")
fcncKinReco_RelativeStatHF = lambda : FCNCKinematicReco(jetVar="RelativeStatHF")
fcncKinReco_PileUpDataMC = lambda : FCNCKinematicReco(jetVar="PileUpDataMC")
fcncKinReco_PileUpPtRef = lambda : FCNCKinematicReco(jetVar="PileUpPtRef")
fcncKinReco_PileUpPtBB = lambda : FCNCKinematicReco(jetVar="PileUpPtBB")
fcncKinReco_PileUpPtEC1 = lambda : FCNCKinematicReco(jetVar="PileUpPtEC1")
fcncKinReco_PileUpPtEC2 = lambda : FCNCKinematicReco(jetVar="PileUpPtEC2")
fcncKinReco_PileUpPtHF = lambda : FCNCKinematicReco(jetVar="PileUpPtHF")
