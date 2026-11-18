import ROOT
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True

import yaml
import os
from ROOT import TLorentzVector
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class FCNCMVAinput(Module, object):
    def __init__(self, *args, **kwargs):
        pass

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        # Branch definitions
        self.out.branch("xsecNorm", "F")
        self.out.branch("MVAinput_Status", "i")

        self.out.branch("MVAinput_WLZL1_dR", "F")
        self.out.branch("MVAinput_WLZL1_dPhi", "F")
        self.out.branch("MVAinput_WLZL1_dEta", "F")

        self.out.branch("MVAinput_WLZL2_dR", "F")
        self.out.branch("MVAinput_WLZL2_dPhi", "F")
        self.out.branch("MVAinput_WLZL2_dEta", "F")

        self.out.branch("MVAinput_ZL1ZL2_dR", "F")
        self.out.branch("MVAinput_ZL1ZL2_dPhi", "F")
        self.out.branch("MVAinput_ZL1ZL2_dEta", "F")

        self.out.branch("MVAinput_Wlepchapt", "F")
        self.out.branch("MVAinput_Wlepchaeta", "F")

        self.out.branch("MVAinput_WJ_DeepJetB", "F")
        self.out.branch("MVAinput_WJ_DeepJetCvL", "F")
        self.out.branch("MVAinput_WJ_DeepJetCvB", "F")
        self.out.branch("MVAinput_WJ_DeepJetQG", "F")
        self.out.branch("MVAinput_WJ_pt", "F")

        self.out.branch("MVAinput_WLWJ_dR", "F")
        self.out.branch("MVAinput_WLWJ_dPhi", "F")
        self.out.branch("MVAinput_WLWJ_dEta", "F")

        self.out.branch("MVAinput_ZL1WJ_dR", "F")
        self.out.branch("MVAinput_ZL1WJ_dPhi", "F")
        self.out.branch("MVAinput_ZL1WJ_dEta", "F")

        self.out.branch("MVAinput_ZL2WJ_dR", "F")
        self.out.branch("MVAinput_ZL2WJ_dPhi", "F")
        self.out.branch("MVAinput_ZL2WJ_dEta", "F")

        self.out.branch("MVAinput_Z_WJ_dR", "F")
        self.out.branch("MVAinput_Z_WJ_dPhi", "F")
        self.out.branch("MVAinput_Z_WJ_dEta", "F")

        self.out.branch("MVAinput_ZJ_DeepJetB", "F")
        self.out.branch("MVAinput_ZJ_DeepJetCvL", "F")
        self.out.branch("MVAinput_ZJ_DeepJetCvB", "F")
        self.out.branch("MVAinput_ZJ_DeepJetQG", "F")
        self.out.branch("MVAinput_ZJ_pt", "F")

        self.out.branch("MVAinput_WLZJ_dR", "F")
        self.out.branch("MVAinput_WLZJ_dPhi", "F")
        self.out.branch("MVAinput_WLZJ_dEta", "F")

        self.out.branch("MVAinput_ZL1ZJ_dR", "F")
        self.out.branch("MVAinput_ZL1ZJ_dPhi", "F")
        self.out.branch("MVAinput_ZL1ZJ_dEta", "F")

        self.out.branch("MVAinput_ZL2ZJ_dR", "F")
        self.out.branch("MVAinput_ZL2ZJ_dPhi", "F")
        self.out.branch("MVAinput_ZL2ZJ_dEta", "F")

        self.out.branch("MVAinput_Z_ZJ_dR", "F")
        self.out.branch("MVAinput_Z_ZJ_dPhi", "F")
        self.out.branch("MVAinput_Z_ZJ_dEta", "F")

        self.out.branch("MVAinput_WJZJ_dR", "F")
        self.out.branch("MVAinput_WJZJ_dPhi", "F")
        self.out.branch("MVAinput_WJZJ_dEta", "F")

        # ---------------------------------------------------------
        # OPTIMIZATION: Load YAML and calculate xsecNorm ONCE per file
        # ---------------------------------------------------------
        info = {}
        cmssw_base = os.environ.get("CMSSW_BASE", "")
        xsecfile = "%s/src/TZWi/TopAnalysis/test/fcncTriLepton/config/crosssection.yaml" % cmssw_base
        
        info.update(yaml.safe_load(open(xsecfile)))
        info.update(yaml.safe_load(open("%s/src/TZWi/TopAnalysis/test/fcncTriLepton/NanoAOD/2018/MC.RunIISummer18.central.yaml" % cmssw_base)))
        info['dataset'].update(yaml.safe_load(open("%s/src/TZWi/TopAnalysis/test/fcncTriLepton/NanoAOD/2018/MC.RunIISummer18.rareprocess.yaml" % cmssw_base))['dataset'])

        rootpath = ROOT.gDirectory.GetPath()
        try:
            rootsample1 = rootpath.split('/')[3]
            rootsample = rootsample1.split('.')[0]
        except IndexError:
            rootsample = "Unknown"

        xsecitem = list(info['crosssection'].items())
        entryitem = list(info['Entries'].items())
        dataitem = list(info['dataset'].items())

        fileset = []
        fileroute = []

        for i, items in enumerate(dataitem):
            fileset0 = items[0]
            fileset1 = fileset0.split('.')[1]
            fileset.append(fileset1)
            fileroutekey0 = list(items[1].keys())
            fileroute.append(fileroutekey0)

        findsample = ""
        for j, routes in enumerate(fileroute):
            for k, routeitem in enumerate(routes):
                if rootsample in routeitem:
                    findsample = fileset[j]

        self.targetxsecweight = 1.0
        targetentry = 1
        
        # Check if the sample is Data
        is_data = rootsample in ["EGamma", "DoubleMuon", "SingleMuon", "MuonEG", "Tau"]
        
        if not is_data and findsample:
            targetxsec = 1.0
            for i, xsecslot in enumerate(xsecitem):
                if findsample == xsecslot[0]:
                    targetxsec = xsecslot[1]
            for i, entryslot in enumerate(entryitem):
                if findsample == entryslot[0]:
                    targetentry = entryslot[1]

            self.targetxsecweight = (float(targetxsec) / float(targetentry)) * 59800.0

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        # Use the pre-calculated cross-section weight
        self.out.fillBranch("xsecNorm", self.targetxsecweight)

        # Variables initialize
        Status      = 0

        WLZL1_dR    = 0
        WLZL1_dPhi  = 0
        WLZL1_dEta  = 0

        WLZL2_dR    = 0
        WLZL2_dPhi  = 0
        WLZL2_dEta  = 0

        ZL1ZL2_dR   = 0
        ZL1ZL2_dPhi = 0
        ZL1ZL2_dEta = 0

        Wlepchapt   = 0
        Wlepchaeta  = 0

        WJ_DeepJetB   = 0
        WJ_DeepJetCvL = 0
        WJ_DeepJetCvB = 0
        WJ_DeepJetQG  = 0
        WJ_pt       = 0

        WLWJ_dR     = 0
        WLWJ_dPhi   = 0
        WLWJ_dEta   = 0

        ZL1WJ_dR    = 0
        ZL1WJ_dPhi  = 0
        ZL1WJ_dEta  = 0

        ZL2WJ_dR    = 0
        ZL2WJ_dPhi  = 0
        ZL2WJ_dEta  = 0

        Z_WJ_dR     = 0
        Z_WJ_dPhi   = 0
        Z_WJ_dEta   = 0

        ZJ_DeepJetB   = 0
        ZJ_DeepJetCvL = 0
        ZJ_DeepJetCvB = 0
        ZJ_DeepJetQG  = 0
        ZJ_pt       = 0

        WLZJ_dR     = 0
        WLZJ_dPhi   = 0
        WLZJ_dEta   = 0

        ZL1ZJ_dR    = 0
        ZL1ZJ_dPhi  = 0
        ZL1ZJ_dEta  = 0

        ZL2ZJ_dR    = 0
        ZL2ZJ_dPhi  = 0
        ZL2ZJ_dEta  = 0

        Z_ZJ_dR     = 0
        Z_ZJ_dPhi   = 0
        Z_ZJ_dEta   = 0

        WJZJ_dR     = 0
        WJZJ_dPhi   = 0
        WJZJ_dEta   = 0

        # CRITICAL FIX: Safe attribute checks. 
        # Verify FCNCKinematicReco successfully reconstructed the objects before using them.
        has_Wlep = hasattr(event._tree, "b_out_Wlepvar")
        has_Zlep1 = hasattr(event._tree, "b_out_Zlep1var")
        has_Zlep2 = hasattr(event._tree, "b_out_Zlep2var")
        has_leptons = has_Wlep and has_Zlep1 and has_Zlep2

        has_W_jet = hasattr(event._tree, "b_out_jet_W")
        has_Z_jet = hasattr(event._tree, "b_out_jet_Z")

        # Basic event selection for WZCR/TTCR
        if getattr(event._tree, "b_out_GoodLeptonCode", 0) == 111 and has_leptons:

            if event._tree.b_out_Lepton1_pdgId > 0:
                Wlepcha = -1
            else:
                Wlepcha = +1

            # Lepton vs Lepton
            WLZL1_dR    = abs(event._tree.b_out_Wlepvar.DeltaR(event._tree.b_out_Zlep1var))
            WLZL1_dPhi  = abs(event._tree.b_out_Wlepvar.DeltaPhi(event._tree.b_out_Zlep1var))
            WLZL1_dEta  = abs(event._tree.b_out_Wlepvar.Eta() - event._tree.b_out_Zlep1var.Eta())

            WLZL2_dR    = abs(event._tree.b_out_Wlepvar.DeltaR(event._tree.b_out_Zlep2var))
            WLZL2_dPhi  = abs(event._tree.b_out_Wlepvar.DeltaPhi(event._tree.b_out_Zlep2var))
            WLZL2_dEta  = abs(event._tree.b_out_Wlepvar.Eta() - event._tree.b_out_Zlep2var.Eta())

            ZL1ZL2_dR   = abs(event._tree.b_out_Zlep1var.DeltaR(event._tree.b_out_Zlep2var))
            ZL1ZL2_dPhi = abs(event._tree.b_out_Zlep1var.DeltaPhi(event._tree.b_out_Zlep2var))
            ZL1ZL2_dEta = abs(event._tree.b_out_Zlep1var.Eta() - event._tree.b_out_Zlep2var.Eta())

            Wlepchapt   = Wlepcha * event._tree.b_out_Wlepvar.Pt()
            Wlepchaeta  = Wlepcha * abs(event._tree.b_out_Wlepvar.Eta())

        if getattr(event._tree, "b_out_GoodLeptonCode", 0) == 111 and getattr(event._tree, "b_out_nGoodJet", 0) >= 1 and getattr(event._tree, "b_out_nBjet", 0) >= 1 and has_leptons and has_W_jet:

            # W Jet & Lepton vs W Jet
            # Safe to read DeepFlav attributes here because we checked has_W_jet
            WJ_DeepJetB   = getattr(event._tree, "b_out_jet_W_DeepFlavB", 0)
            WJ_DeepJetCvL = getattr(event._tree, "b_out_jet_W_DeepFlavCvL", 0)
            WJ_DeepJetCvB = getattr(event._tree, "b_out_jet_W_DeepFlavCvB", 0)
            WJ_DeepJetQG  = getattr(event._tree, "b_out_jet_W_DeepFlavQG", 0)
            WJ_pt       = event._tree.b_out_jet_W.Pt()

            WLWJ_dR     = abs(event._tree.b_out_Wlepvar.DeltaR(event._tree.b_out_jet_W))
            WLWJ_dPhi   = abs(event._tree.b_out_Wlepvar.DeltaPhi(event._tree.b_out_jet_W))
            WLWJ_dEta   = abs(event._tree.b_out_Wlepvar.Eta() - event._tree.b_out_jet_W.Eta())

            ZL1WJ_dR    = abs(event._tree.b_out_Zlep1var.DeltaR(event._tree.b_out_jet_W))
            ZL1WJ_dPhi  = abs(event._tree.b_out_Zlep1var.DeltaPhi(event._tree.b_out_jet_W))
            ZL1WJ_dEta  = abs(event._tree.b_out_Zlep1var.Eta() - event._tree.b_out_jet_W.Eta())

            ZL2WJ_dR    = abs(event._tree.b_out_Zlep2var.DeltaR(event._tree.b_out_jet_W))
            ZL2WJ_dPhi  = abs(event._tree.b_out_Zlep2var.DeltaPhi(event._tree.b_out_jet_W))
            ZL2WJ_dEta  = abs(event._tree.b_out_Zlep2var.Eta() - event._tree.b_out_jet_W.Eta())

            Z_WJ_dR     = abs((event._tree.b_out_Zlep1var+event._tree.b_out_Zlep2var).DeltaR(event._tree.b_out_jet_W))
            Z_WJ_dPhi   = abs((event._tree.b_out_Zlep1var+event._tree.b_out_Zlep2var).DeltaPhi(event._tree.b_out_jet_W))
            Z_WJ_dEta   = abs((event._tree.b_out_Zlep1var+event._tree.b_out_Zlep2var).Eta() - event._tree.b_out_jet_W.Eta())

        if getattr(event._tree, "b_out_GoodLeptonCode", 0) == 111 and getattr(event._tree, "b_out_nGoodJet", 0) >= 2 and getattr(event._tree, "b_out_nBjet", 0) >= 1 and has_leptons and has_Z_jet:

            # Z Jet & Lepton vs Z Jet
            ZJ_DeepJetB   = getattr(event._tree, "b_out_jet_Z_DeepFlavB", 0)
            ZJ_DeepJetCvL = getattr(event._tree, "b_out_jet_Z_DeepFlavCvL", 0)
            ZJ_DeepJetCvB = getattr(event._tree, "b_out_jet_Z_DeepFlavCvB", 0)
            ZJ_DeepJetQG  = getattr(event._tree, "b_out_jet_Z_DeepFlavQG", 0)
            ZJ_pt       = event._tree.b_out_jet_Z.Pt()

            WLZJ_dR     = abs(event._tree.b_out_Wlepvar.DeltaR(event._tree.b_out_jet_Z))
            WLZJ_dPhi   = abs(event._tree.b_out_Wlepvar.DeltaPhi(event._tree.b_out_jet_Z))
            WLZJ_dEta   = abs(event._tree.b_out_Wlepvar.Eta() - event._tree.b_out_jet_Z.Eta())

            ZL1ZJ_dR    = abs(event._tree.b_out_Zlep1var.DeltaR(event._tree.b_out_jet_Z))
            ZL1ZJ_dPhi  = abs(event._tree.b_out_Zlep1var.DeltaPhi(event._tree.b_out_jet_Z))
            ZL1ZJ_dEta  = abs(event._tree.b_out_Zlep1var.Eta() - event._tree.b_out_jet_Z.Eta())

            ZL2ZJ_dR    = abs(event._tree.b_out_Zlep2var.DeltaR(event._tree.b_out_jet_Z))
            ZL2ZJ_dPhi  = abs(event._tree.b_out_Zlep2var.DeltaPhi(event._tree.b_out_jet_Z))
            ZL2ZJ_dEta  = abs(event._tree.b_out_Zlep2var.Eta() - event._tree.b_out_jet_Z.Eta())

            Z_ZJ_dR     = abs((event._tree.b_out_Zlep1var+event._tree.b_out_Zlep2var).DeltaR(event._tree.b_out_jet_Z))
            Z_ZJ_dPhi   = abs((event._tree.b_out_Zlep1var+event._tree.b_out_Zlep2var).DeltaPhi(event._tree.b_out_jet_Z))
            Z_ZJ_dEta   = abs((event._tree.b_out_Zlep1var+event._tree.b_out_Zlep2var).Eta() - event._tree.b_out_jet_Z.Eta())

            # W Jet vs Z Jet
            # For this we need BOTH has_W_jet and has_Z_jet
            if has_W_jet:
                WJZJ_dR     = abs(event._tree.b_out_jet_W.DeltaR(event._tree.b_out_jet_Z))
                WJZJ_dPhi   = abs(event._tree.b_out_jet_W.DeltaPhi(event._tree.b_out_jet_Z))
                WJZJ_dEta   = abs(event._tree.b_out_jet_W.Eta() - event._tree.b_out_jet_Z.Eta())

        self.out.fillBranch("MVAinput_Status", Status)

        self.out.fillBranch("MVAinput_WLZL1_dR", WLZL1_dR)
        self.out.fillBranch("MVAinput_WLZL1_dPhi", WLZL1_dPhi)
        self.out.fillBranch("MVAinput_WLZL1_dEta", WLZL1_dEta)

        self.out.fillBranch("MVAinput_WLZL2_dR", WLZL2_dR)
        self.out.fillBranch("MVAinput_WLZL2_dPhi", WLZL2_dPhi)
        self.out.fillBranch("MVAinput_WLZL2_dEta", WLZL2_dEta)

        self.out.fillBranch("MVAinput_ZL1ZL2_dR", ZL1ZL2_dR)
        self.out.fillBranch("MVAinput_ZL1ZL2_dPhi", ZL1ZL2_dPhi)
        self.out.fillBranch("MVAinput_ZL1ZL2_dEta", ZL1ZL2_dEta)

        self.out.fillBranch("MVAinput_Wlepchapt", Wlepchapt)
        self.out.fillBranch("MVAinput_Wlepchaeta", Wlepchaeta)

        self.out.fillBranch("MVAinput_WJ_DeepJetB", WJ_DeepJetB)
        self.out.fillBranch("MVAinput_WJ_DeepJetCvL", WJ_DeepJetCvL)
        self.out.fillBranch("MVAinput_WJ_DeepJetCvB", WJ_DeepJetCvB)
        self.out.fillBranch("MVAinput_WJ_DeepJetQG", WJ_DeepJetQG)
        self.out.fillBranch("MVAinput_WJ_pt", WJ_pt)

        self.out.fillBranch("MVAinput_WLWJ_dR", WLWJ_dR)
        self.out.fillBranch("MVAinput_WLWJ_dPhi", WLWJ_dPhi)
        self.out.fillBranch("MVAinput_WLWJ_dEta", WLWJ_dEta)

        self.out.fillBranch("MVAinput_ZL1WJ_dR", ZL1WJ_dR)
        self.out.fillBranch("MVAinput_ZL1WJ_dPhi", ZL1WJ_dPhi)
        self.out.fillBranch("MVAinput_ZL1WJ_dEta", ZL1WJ_dEta)

        self.out.fillBranch("MVAinput_ZL2WJ_dR", ZL2WJ_dR)
        self.out.fillBranch("MVAinput_ZL2WJ_dPhi", ZL2WJ_dPhi)
        self.out.fillBranch("MVAinput_ZL2WJ_dEta", ZL2WJ_dEta)

        self.out.fillBranch("MVAinput_Z_WJ_dR", Z_WJ_dR)
        self.out.fillBranch("MVAinput_Z_WJ_dPhi", Z_WJ_dPhi)
        self.out.fillBranch("MVAinput_Z_WJ_dEta", Z_WJ_dEta)

        self.out.fillBranch("MVAinput_ZJ_DeepJetB", ZJ_DeepJetB)
        self.out.fillBranch("MVAinput_ZJ_DeepJetCvL", ZJ_DeepJetCvL)
        self.out.fillBranch("MVAinput_ZJ_DeepJetCvB", ZJ_DeepJetCvB)
        self.out.fillBranch("MVAinput_ZJ_DeepJetQG", ZJ_DeepJetQG)
        self.out.fillBranch("MVAinput_ZJ_pt", ZJ_pt)

        self.out.fillBranch("MVAinput_WLZJ_dR", WLZJ_dR)
        self.out.fillBranch("MVAinput_WLZJ_dPhi", WLZJ_dPhi)
        self.out.fillBranch("MVAinput_WLZJ_dEta", WLZJ_dEta)

        self.out.fillBranch("MVAinput_ZL1ZJ_dR", ZL1ZJ_dR)
        self.out.fillBranch("MVAinput_ZL1ZJ_dPhi", ZL1ZJ_dPhi)
        self.out.fillBranch("MVAinput_ZL1ZJ_dEta", ZL1ZJ_dEta)

        self.out.fillBranch("MVAinput_ZL2ZJ_dR", ZL2ZJ_dR)
        self.out.fillBranch("MVAinput_ZL2ZJ_dPhi", ZL2ZJ_dPhi)
        self.out.fillBranch("MVAinput_ZL2ZJ_dEta", ZL2ZJ_dEta)

        self.out.fillBranch("MVAinput_Z_ZJ_dR", Z_ZJ_dR)
        self.out.fillBranch("MVAinput_Z_ZJ_dPhi", Z_ZJ_dPhi)
        self.out.fillBranch("MVAinput_Z_ZJ_dEta", Z_ZJ_dEta)

        self.out.fillBranch("MVAinput_WJZJ_dR", WJZJ_dR)
        self.out.fillBranch("MVAinput_WJZJ_dPhi", WJZJ_dPhi)
        self.out.fillBranch("MVAinput_WJZJ_dEta", WJZJ_dEta)

        return True

fcncMVAinput = lambda: FCNCMVAinput()
