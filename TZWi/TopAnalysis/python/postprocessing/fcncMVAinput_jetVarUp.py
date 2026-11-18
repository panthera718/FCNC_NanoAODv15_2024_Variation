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
        self.jetVar = kwargs.get("jetVar")
        pass

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        jetVar = self.jetVar

        # Branch definitions
        self.out.branch("xsecNorm_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_Status_%sUp" % jetVar, "i")

        self.out.branch("MVAinput_WLZL1_dR_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_WLZL1_dPhi_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_WLZL1_dEta_%sUp" % jetVar, "F")

        self.out.branch("MVAinput_WLZL2_dR_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_WLZL2_dPhi_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_WLZL2_dEta_%sUp" % jetVar, "F")

        self.out.branch("MVAinput_ZL1ZL2_dR_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_ZL1ZL2_dPhi_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_ZL1ZL2_dEta_%sUp" % jetVar, "F")

        self.out.branch("MVAinput_Wlepchapt_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_Wlepchaeta_%sUp" % jetVar, "F")

        self.out.branch("MVAinput_WJ_DeepJetB_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_WJ_DeepJetCvL_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_WJ_DeepJetCvB_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_WJ_DeepJetQG_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_WJ_pt_%sUp" % jetVar, "F")

        self.out.branch("MVAinput_WLWJ_dR_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_WLWJ_dPhi_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_WLWJ_dEta_%sUp" % jetVar, "F")

        self.out.branch("MVAinput_ZL1WJ_dR_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_ZL1WJ_dPhi_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_ZL1WJ_dEta_%sUp" % jetVar, "F")

        self.out.branch("MVAinput_ZL2WJ_dR_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_ZL2WJ_dPhi_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_ZL2WJ_dEta_%sUp" % jetVar, "F")

        self.out.branch("MVAinput_Z_WJ_dR_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_Z_WJ_dPhi_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_Z_WJ_dEta_%sUp" % jetVar, "F")

        self.out.branch("MVAinput_ZJ_DeepJetB_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_ZJ_DeepJetCvL_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_ZJ_DeepJetCvB_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_ZJ_DeepJetQG_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_ZJ_pt_%sUp" % jetVar, "F")

        self.out.branch("MVAinput_WLZJ_dR_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_WLZJ_dPhi_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_WLZJ_dEta_%sUp" % jetVar, "F")

        self.out.branch("MVAinput_ZL1ZJ_dR_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_ZL1ZJ_dPhi_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_ZL1ZJ_dEta_%sUp" % jetVar, "F")

        self.out.branch("MVAinput_ZL2ZJ_dR_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_ZL2ZJ_dPhi_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_ZL2ZJ_dEta_%sUp" % jetVar, "F")

        self.out.branch("MVAinput_Z_ZJ_dR_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_Z_ZJ_dPhi_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_Z_ZJ_dEta_%sUp" % jetVar, "F")

        self.out.branch("MVAinput_WJZJ_dR_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_WJZJ_dPhi_%sUp" % jetVar, "F")
        self.out.branch("MVAinput_WJZJ_dEta_%sUp" % jetVar, "F")

        # ---------------------------------------------------------
        # OPTIMIZATION: Load YAML and calculate xsecNorm ONCE per file
        # ---------------------------------------------------------
        info = {}
        cmssw_base = os.environ.get("CMSSW_BASE", "")
        xsecfile = "%s/src/TZWi/TopAnalysis/test/fcncTriLepton/config/crosssection.yaml" % cmssw_base
        
        info.update(yaml.safe_load(open(xsecfile)))
        info.update(yaml.safe_load(open("%s/src/TZWi/TopAnalysis/test/fcncTriLepton/NanoAOD/2016/MC.RunIISummer16.central.yaml" % cmssw_base)))
        info['dataset'].update(yaml.safe_load(open("%s/src/TZWi/TopAnalysis/test/fcncTriLepton/NanoAOD/2016/MC.RunIISummer16.rareprocess.yaml" % cmssw_base))['dataset'])

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
        is_data = rootsample in ["DoubleEG", "SingleElectron", "DoubleMuon", "SingleMuon", "MuonEG", "Tau"]
        
        if not is_data and findsample:
            targetxsec = 1.0
            for i, xsecslot in enumerate(xsecitem):
                if findsample == xsecslot[0]:
                    targetxsec = xsecslot[1]
            for i, entryslot in enumerate(entryitem):
                if findsample == entryslot[0]:
                    targetentry = entryslot[1]

            self.targetxsecweight = (float(targetxsec) / float(targetentry)) * 19500.0

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        jetVar = self.jetVar

        # Use the pre-calculated cross-section weight
        self.out.fillBranch("xsecNorm_%sUp" % jetVar, self.targetxsecweight)

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
        has_Wlep = hasattr(event._tree, "b_out_Wlepvar_%sUp" % jetVar)
        has_Zlep1 = hasattr(event._tree, "b_out_Zlep1var_%sUp" % jetVar)
        has_Zlep2 = hasattr(event._tree, "b_out_Zlep2var_%sUp" % jetVar)
        has_leptons = has_Wlep and has_Zlep1 and has_Zlep2

        has_W_jet = hasattr(event._tree, "b_out_jet_W_%sUp" % jetVar)
        has_Z_jet = hasattr(event._tree, "b_out_jet_Z_%sUp" % jetVar)

        # Basic event selection for WZCR/TTCR
        if getattr(event._tree, "b_out_GoodLeptonCode", 0) == 111 and has_leptons:

            if event._tree.b_out_Lepton1_pdgId > 0:
                Wlepcha = -1
            else:
                Wlepcha = +1

            b_out_Wlepvar = getattr(event._tree, "b_out_Wlepvar_%sUp" % jetVar)
            b_out_Zlep1var = getattr(event._tree, "b_out_Zlep1var_%sUp" % jetVar)
            b_out_Zlep2var = getattr(event._tree, "b_out_Zlep2var_%sUp" % jetVar)

            # Lepton vs Lepton
            WLZL1_dR    = abs(b_out_Wlepvar.DeltaR(b_out_Zlep1var))
            WLZL1_dPhi  = abs(b_out_Wlepvar.DeltaPhi(b_out_Zlep1var))
            WLZL1_dEta  = abs(b_out_Wlepvar.Eta() - b_out_Zlep1var.Eta())

            WLZL2_dR    = abs(b_out_Wlepvar.DeltaR(b_out_Zlep2var))
            WLZL2_dPhi  = abs(b_out_Wlepvar.DeltaPhi(b_out_Zlep2var))
            WLZL2_dEta  = abs(b_out_Wlepvar.Eta() - b_out_Zlep2var.Eta())

            ZL1ZL2_dR   = abs(b_out_Zlep1var.DeltaR(b_out_Zlep2var))
            ZL1ZL2_dPhi = abs(b_out_Zlep1var.DeltaPhi(b_out_Zlep2var))
            ZL1ZL2_dEta = abs(b_out_Zlep1var.Eta() - b_out_Zlep2var.Eta())

            Wlepchapt   = Wlepcha * b_out_Wlepvar.Pt()
            Wlepchaeta  = Wlepcha * abs(b_out_Wlepvar.Eta())

        if getattr(event._tree, "b_out_GoodLeptonCode", 0) == 111 and getattr(event._tree, "b_out_nGoodJet_%sUp" % jetVar, 0) >= 1 and getattr(event._tree, "b_out_nBjet_%sUp" % jetVar, 0) >= 1 and has_leptons and has_W_jet:

            b_out_Wlepvar = getattr(event._tree, "b_out_Wlepvar_%sUp" % jetVar)
            b_out_Zlep1var = getattr(event._tree, "b_out_Zlep1var_%sUp" % jetVar)
            b_out_Zlep2var = getattr(event._tree, "b_out_Zlep2var_%sUp" % jetVar)
            b_out_jet_W = getattr(event._tree, "b_out_jet_W_%sUp" % jetVar)

            # W Jet & Lepton vs W Jet
            # Safe to read DeepFlav attributes here because we checked has_W_jet
            WJ_DeepJetB   = getattr(event._tree, "b_out_jet_W_DeepFlavB_%sUp" % jetVar, 0)
            WJ_DeepJetCvL = getattr(event._tree, "b_out_jet_W_DeepFlavCvL_%sUp" % jetVar, 0)
            WJ_DeepJetCvB = getattr(event._tree, "b_out_jet_W_DeepFlavCvB_%sUp" % jetVar, 0)
            WJ_DeepJetQG  = getattr(event._tree, "b_out_jet_W_DeepFlavQG_%sUp" % jetVar, 0)
            WJ_pt       = b_out_jet_W.Pt()

            WLWJ_dR     = abs(b_out_Wlepvar.DeltaR(b_out_jet_W))
            WLWJ_dPhi   = abs(b_out_Wlepvar.DeltaPhi(b_out_jet_W))
            WLWJ_dEta   = abs(b_out_Wlepvar.Eta() - b_out_jet_W.Eta())

            ZL1WJ_dR    = abs(b_out_Zlep1var.DeltaR(b_out_jet_W))
            ZL1WJ_dPhi  = abs(b_out_Zlep1var.DeltaPhi(b_out_jet_W))
            ZL1WJ_dEta  = abs(b_out_Zlep1var.Eta() - b_out_jet_W.Eta())

            ZL2WJ_dR    = abs(b_out_Zlep2var.DeltaR(b_out_jet_W))
            ZL2WJ_dPhi  = abs(b_out_Zlep2var.DeltaPhi(b_out_jet_W))
            ZL2WJ_dEta  = abs(b_out_Zlep2var.Eta() - b_out_jet_W.Eta())

            Z_WJ_dR     = abs((b_out_Zlep1var+b_out_Zlep2var).DeltaR(b_out_jet_W))
            Z_WJ_dPhi   = abs((b_out_Zlep1var+b_out_Zlep2var).DeltaPhi(b_out_jet_W))
            Z_WJ_dEta   = abs((b_out_Zlep1var+b_out_Zlep2var).Eta() - b_out_jet_W.Eta())

        if getattr(event._tree, "b_out_GoodLeptonCode", 0) == 111 and getattr(event._tree, "b_out_nGoodJet_%sUp" % jetVar, 0) >= 2 and getattr(event._tree, "b_out_nBjet_%sUp" % jetVar, 0) >= 1 and has_leptons and has_Z_jet:

            b_out_Wlepvar = getattr(event._tree, "b_out_Wlepvar_%sUp" % jetVar)
            b_out_Zlep1var = getattr(event._tree, "b_out_Zlep1var_%sUp" % jetVar)
            b_out_Zlep2var = getattr(event._tree, "b_out_Zlep2var_%sUp" % jetVar)
            b_out_jet_W = getattr(event._tree, "b_out_jet_W_%sUp" % jetVar)
            b_out_jet_Z = getattr(event._tree, "b_out_jet_Z_%sUp" % jetVar)

            # Z Jet & Lepton vs Z Jet
            ZJ_DeepJetB   = getattr(event._tree, "b_out_jet_Z_DeepFlavB_%sUp" % jetVar, 0)
            ZJ_DeepJetCvL = getattr(event._tree, "b_out_jet_Z_DeepFlavCvL_%sUp" % jetVar, 0)
            ZJ_DeepJetCvB = getattr(event._tree, "b_out_jet_Z_DeepFlavCvB_%sUp" % jetVar, 0)
            ZJ_DeepJetQG  = getattr(event._tree, "b_out_jet_Z_DeepFlavQG_%sUp" % jetVar, 0)
            ZJ_pt       = b_out_jet_Z.Pt()

            WLZJ_dR     = abs(b_out_Wlepvar.DeltaR(b_out_jet_Z))
            WLZJ_dPhi   = abs(b_out_Wlepvar.DeltaPhi(b_out_jet_Z))
            WLZJ_dEta   = abs(b_out_Wlepvar.Eta() - b_out_jet_Z.Eta())

            ZL1ZJ_dR    = abs(b_out_Zlep1var.DeltaR(b_out_jet_Z))
            ZL1ZJ_dPhi  = abs(b_out_Zlep1var.DeltaPhi(b_out_jet_Z))
            ZL1ZJ_dEta  = abs(b_out_Zlep1var.Eta() - b_out_jet_Z.Eta())

            ZL2ZJ_dR    = abs(b_out_Zlep2var.DeltaR(b_out_jet_Z))
            ZL2ZJ_dPhi  = abs(b_out_Zlep2var.DeltaPhi(b_out_jet_Z))
            ZL2ZJ_dEta  = abs(b_out_Zlep2var.Eta() - b_out_jet_Z.Eta())

            Z_ZJ_dR     = abs((b_out_Zlep1var+b_out_Zlep2var).DeltaR(b_out_jet_Z))
            Z_ZJ_dPhi   = abs((b_out_Zlep1var+b_out_Zlep2var).DeltaPhi(b_out_jet_Z))
            Z_ZJ_dEta   = abs((b_out_Zlep1var+b_out_Zlep2var).Eta() - b_out_jet_Z.Eta())

            # W Jet vs Z Jet
            # For this we need BOTH has_W_jet and has_Z_jet
            if has_W_jet:
                WJZJ_dR     = abs(b_out_jet_W.DeltaR(b_out_jet_Z))
                WJZJ_dPhi   = abs(b_out_jet_W.DeltaPhi(b_out_jet_Z))
                WJZJ_dEta   = abs(b_out_jet_W.Eta() - b_out_jet_Z.Eta())

        self.out.fillBranch("MVAinput_Status_%sUp" % jetVar, Status)

        self.out.fillBranch("MVAinput_WLZL1_dR_%sUp" % jetVar, WLZL1_dR)
        self.out.fillBranch("MVAinput_WLZL1_dPhi_%sUp" % jetVar, WLZL1_dPhi)
        self.out.fillBranch("MVAinput_WLZL1_dEta_%sUp" % jetVar, WLZL1_dEta)

        self.out.fillBranch("MVAinput_WLZL2_dR_%sUp" % jetVar, WLZL2_dR)
        self.out.fillBranch("MVAinput_WLZL2_dPhi_%sUp" % jetVar, WLZL2_dPhi)
        self.out.fillBranch("MVAinput_WLZL2_dEta_%sUp" % jetVar, WLZL2_dEta)

        self.out.fillBranch("MVAinput_ZL1ZL2_dR_%sUp" % jetVar, ZL1ZL2_dR)
        self.out.fillBranch("MVAinput_ZL1ZL2_dPhi_%sUp" % jetVar, ZL1ZL2_dPhi)
        self.out.fillBranch("MVAinput_ZL1ZL2_dEta_%sUp" % jetVar, ZL1ZL2_dEta)

        self.out.fillBranch("MVAinput_Wlepchapt_%sUp" % jetVar, Wlepchapt)
        self.out.fillBranch("MVAinput_Wlepchaeta_%sUp" % jetVar, Wlepchaeta)

        self.out.fillBranch("MVAinput_WJ_DeepJetB_%sUp" % jetVar, WJ_DeepJetB)
        self.out.fillBranch("MVAinput_WJ_DeepJetCvL_%sUp" % jetVar, WJ_DeepJetCvL)
        self.out.fillBranch("MVAinput_WJ_DeepJetCvB_%sUp" % jetVar, WJ_DeepJetCvB)
        self.out.fillBranch("MVAinput_WJ_DeepJetQG_%sUp" % jetVar, WJ_DeepJetQG)
        self.out.fillBranch("MVAinput_WJ_pt_%sUp" % jetVar, WJ_pt)

        self.out.fillBranch("MVAinput_WLWJ_dR_%sUp" % jetVar, WLWJ_dR)
        self.out.fillBranch("MVAinput_WLWJ_dPhi_%sUp" % jetVar, WLWJ_dPhi)
        self.out.fillBranch("MVAinput_WLWJ_dEta_%sUp" % jetVar, WLWJ_dEta)

        self.out.fillBranch("MVAinput_ZL1WJ_dR_%sUp" % jetVar, ZL1WJ_dR)
        self.out.fillBranch("MVAinput_ZL1WJ_dPhi_%sUp" % jetVar, ZL1WJ_dPhi)
        self.out.fillBranch("MVAinput_ZL1WJ_dEta_%sUp" % jetVar, ZL1WJ_dEta)

        self.out.fillBranch("MVAinput_ZL2WJ_dR_%sUp" % jetVar, ZL2WJ_dR)
        self.out.fillBranch("MVAinput_ZL2WJ_dPhi_%sUp" % jetVar, ZL2WJ_dPhi)
        self.out.fillBranch("MVAinput_ZL2WJ_dEta_%sUp" % jetVar, ZL2WJ_dEta)

        self.out.fillBranch("MVAinput_Z_WJ_dR_%sUp" % jetVar, Z_WJ_dR)
        self.out.fillBranch("MVAinput_Z_WJ_dPhi_%sUp" % jetVar, Z_WJ_dPhi)
        self.out.fillBranch("MVAinput_Z_WJ_dEta_%sUp" % jetVar, Z_WJ_dEta)

        self.out.fillBranch("MVAinput_ZJ_DeepJetB_%sUp" % jetVar, ZJ_DeepJetB)
        self.out.fillBranch("MVAinput_ZJ_DeepJetCvL_%sUp" % jetVar, ZJ_DeepJetCvL)
        self.out.fillBranch("MVAinput_ZJ_DeepJetCvB_%sUp" % jetVar, ZJ_DeepJetCvB)
        self.out.fillBranch("MVAinput_ZJ_DeepJetQG_%sUp" % jetVar, ZJ_DeepJetQG)
        self.out.fillBranch("MVAinput_ZJ_pt_%sUp" % jetVar, ZJ_pt)

        self.out.fillBranch("MVAinput_WLZJ_dR_%sUp" % jetVar, WLZJ_dR)
        self.out.fillBranch("MVAinput_WLZJ_dPhi_%sUp" % jetVar, WLZJ_dPhi)
        self.out.fillBranch("MVAinput_WLZJ_dEta_%sUp" % jetVar, WLZJ_dEta)

        self.out.fillBranch("MVAinput_ZL1ZJ_dR_%sUp" % jetVar, ZL1ZJ_dR)
        self.out.fillBranch("MVAinput_ZL1ZJ_dPhi_%sUp" % jetVar, ZL1ZJ_dPhi)
        self.out.fillBranch("MVAinput_ZL1ZJ_dEta_%sUp" % jetVar, ZL1ZJ_dEta)

        self.out.fillBranch("MVAinput_ZL2ZJ_dR_%sUp" % jetVar, ZL2ZJ_dR)
        self.out.fillBranch("MVAinput_ZL2ZJ_dPhi_%sUp" % jetVar, ZL2ZJ_dPhi)
        self.out.fillBranch("MVAinput_ZL2ZJ_dEta_%sUp" % jetVar, ZL2ZJ_dEta)

        self.out.fillBranch("MVAinput_Z_ZJ_dR_%sUp" % jetVar, Z_ZJ_dR)
        self.out.fillBranch("MVAinput_Z_ZJ_dPhi_%sUp" % jetVar, Z_ZJ_dPhi)
        self.out.fillBranch("MVAinput_Z_ZJ_dEta_%sUp" % jetVar, Z_ZJ_dEta)

        self.out.fillBranch("MVAinput_WJZJ_dR_%sUp" % jetVar, WJZJ_dR)
        self.out.fillBranch("MVAinput_WJZJ_dPhi_%sUp" % jetVar, WJZJ_dPhi)
        self.out.fillBranch("MVAinput_WJZJ_dEta_%sUp" % jetVar, WJZJ_dEta)

        return True


fcncMVAinput_AbsoluteStat = lambda : FCNCMVAinput(jetVar="AbsoluteStat")
fcncMVAinput_AbsoluteScale = lambda : FCNCMVAinput(jetVar="AbsoluteScale")
fcncMVAinput_AbsoluteMPFBias = lambda : FCNCMVAinput(jetVar="AbsoluteMPFBias")
fcncMVAinput_Fragmentation = lambda : FCNCMVAinput(jetVar="Fragmentation")
fcncMVAinput_SinglePionECAL = lambda : FCNCMVAinput(jetVar="SinglePionECAL")
fcncMVAinput_SinglePionHCAL = lambda : FCNCMVAinput(jetVar="SinglePionHCAL")
fcncMVAinput_FlavorQCD = lambda : FCNCMVAinput(jetVar="FlavorQCD")
fcncMVAinput_FlavorPureQuark = lambda : FCNCMVAinput(jetVar="FlavorPureQuark")
fcncMVAinput_FlavorPureCharm = lambda : FCNCMVAinput(jetVar="FlavorPureCharm")
fcncMVAinput_FlavorPureBottom = lambda : FCNCMVAinput(jetVar="FlavorPureBottom")
fcncMVAinput_TimePtEta = lambda : FCNCMVAinput(jetVar="TimePtEta")
fcncMVAinput_RelativeJEREC1 = lambda : FCNCMVAinput(jetVar="RelativeJEREC1")
fcncMVAinput_RelativeJEREC2 = lambda : FCNCMVAinput(jetVar="RelativeJEREC2")
fcncMVAinput_RelativeJERHF = lambda : FCNCMVAinput(jetVar="RelativeJERHF")
fcncMVAinput_RelativePtBB = lambda : FCNCMVAinput(jetVar="RelativePtBB")
fcncMVAinput_RelativePtEC1 = lambda : FCNCMVAinput(jetVar="RelativePtEC1")
fcncMVAinput_RelativePtEC2 = lambda : FCNCMVAinput(jetVar="RelativePtEC2")
fcncMVAinput_RelativePtHF = lambda : FCNCMVAinput(jetVar="RelativePtHF")
fcncMVAinput_RelativeBal = lambda : FCNCMVAinput(jetVar="RelativeBal")
fcncMVAinput_RelativeSample = lambda : FCNCMVAinput(jetVar="RelativeSample")
fcncMVAinput_RelativeFSR = lambda : FCNCMVAinput(jetVar="RelativeFSR")
fcncMVAinput_RelativeStatFSR = lambda : FCNCMVAinput(jetVar="RelativeStatFSR")
fcncMVAinput_RelativeStatEC = lambda : FCNCMVAinput(jetVar="RelativeStatEC")
fcncMVAinput_RelativeStatHF = lambda : FCNCMVAinput(jetVar="RelativeStatHF")
fcncMVAinput_PileUpDataMC = lambda : FCNCMVAinput(jetVar="PileUpDataMC")
fcncMVAinput_PileUpPtRef = lambda : FCNCMVAinput(jetVar="PileUpPtRef")
fcncMVAinput_PileUpPtBB = lambda : FCNCMVAinput(jetVar="PileUpPtBB")
fcncMVAinput_PileUpPtEC1 = lambda : FCNCMVAinput(jetVar="PileUpPtEC1")
fcncMVAinput_PileUpPtEC2 = lambda : FCNCMVAinput(jetVar="PileUpPtEC2")
fcncMVAinput_PileUpPtHF = lambda : FCNCMVAinput(jetVar="PileUpPtHF")
