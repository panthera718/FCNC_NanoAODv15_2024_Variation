import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import os
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class FCNCTriLepton(Module, object):
    def __init__(self, *args, **kwargs):
        self.mode = kwargs.get("mode")
        self.jetVar = kwargs.get("jetVar")

        if "/JetUncertUp_FCNCTriLeptonCppWorker_cc.so" not in  ROOT.gSystem.GetLibraries():
            print("Load C++ FCNCTriLepton worker module")
            base = os.getenv("NANOAODTOOLS_BASE")
            if base:
                ROOT.gROOT.ProcessLine(".L %s/src/JetUncertUp_FCNCTriLeptonCppWorker.cc+O" % base)
            else:
                base = "%s/src/TZWi/TopAnalysis"%os.getenv("CMSSW_BASE")
                ROOT.gSystem.Load("libTZWiTopAnalysis.so")
                ROOT.gROOT.ProcessLine(".L %s/interface/JetUncertUp_FCNCTriLeptonCppWorker.h" % base)
        pass

    def beginJob(self):
        self.worker = ROOT.JetUncertUp_FCNCTriLeptonCppWorker(self.mode)
        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        # 1. Clean up TTreeReader objects to prevent cppyy double-free segfaults
        for attr in list(self.__dict__.keys()):
            if attr.startswith("b_"):
                delattr(self, attr)
                
        # 2. Clean up dynamically attached variables on the PyROOT TTree proxy
        if inputTree and hasattr(inputTree, 'tree') and inputTree.tree():
            tree = inputTree.tree()
            for attr in list(tree.__dict__.keys()):
                if attr.startswith("b_out_"):
                    delattr(tree, attr)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        jetVar = self.jetVar

        self.out.branch("MET_pt_%sUp" % jetVar, "F")
        self.out.branch("MET_phi_%sUp" % jetVar, "F")

        self.out.branch("GoodJet_index_%sUp" % jetVar, "i", lenVar="nGoodJet_%sUp" % jetVar)
        for varName in ["pt", "eta", "phi", "mass", "DeepFlavB"]:
            self.out.branch("GoodJet_%s_%sUp" % (varName, jetVar), "F", lenVar="nGoodJet_%sUp" % jetVar)
        for varName in ["hadronFlavour"]:
            self.out.branch("GoodJet_%s_%sUp" % (varName, jetVar), "I", lenVar="nGoodJet_%sUp" % jetVar)
        for objName in ["ptsortGoodJet1", "ptsortGoodJet2", "ptsortGoodJet3"]:
            for varName in ["pt", "eta", "phi"]:
                self.out.branch("%s_%s_%sUp" % (objName, varName, jetVar), "F")
        
        self.out.branch("nBjet_%sUp" % jetVar, "i")
        self.out.branch("W_MT_%sUp" % jetVar, "F")

        # 2. In the beginFile method, define the new branches with JetVar suffixes:
        self.out.branch("nCjet_%sUp" % jetVar, "i")
        self.out.branch("nUDSjet_%sUp" % jetVar, "i")
        self.out.branch("GoodJet_DeepFlavCvL_%sUp" % jetVar, "F", lenVar="nGoodJet_%sUp" % jetVar)
        self.out.branch("GoodJet_DeepFlavCvB_%sUp" % jetVar, "F", lenVar="nGoodJet_%sUp" % jetVar)
        self.out.branch("GoodJet_DeepFlavQG_%sUp" % jetVar, "F", lenVar="nGoodJet_%sUp" % jetVar)

        self.initReaders(inputTree)
        pass

    def initReaders(self,tree):
        self.b_MET_pt = tree.valueReader("MET_T1_pt")
        self.b_MET_phi = tree.valueReader("MET_T1_phi")

        objName = "Electron"
        setattr(self, "b_n%s" % objName, tree.valueReader("n%s" % objName))
        for varName in ["pt", "eta", "isPFcand", "lostHits", "dxy", "dz", "sip3d",
                        "miniPFRelIso_all", "sieie", "hoe", "eInvMinusPInv", "convVeto",
                        "mvaTOP", "phi", "mass", "charge",
                        "jetPt", "jetEta", "jetPhi", "jetMass",
                        "jetPtRatio", "jetRelIso", "mvaFall17V2noIso_WPL", "pfRelIso03_all", "jetBTagDeepFlavor", "deltaEtaSC"]:
            setattr(self, "b_%s_%s" % (objName, varName), tree.arrayReader("%s_%s" % (objName, varName)))

        objName = "Muon"
        setattr(self, "b_n%s" % objName, tree.valueReader("n%s" % objName))
        for varName in ["pt", "eta", "isPFcand", "isTracker", "isGlobal", "dxy", "dz", "sip3d",
                        "miniPFRelIso_all", "mediumId",
                        "mvaTOP", "phi", "mass", "charge",
                        "jetPt", "jetEta", "jetPhi", "jetMass", "jetPtRatio", "pfRelIso04_all", "jetBTagDeepFlavor"]:
            setattr(self, "b_%s_%s" % (objName, varName), tree.arrayReader("%s_%s" % (objName, varName)))

        objName = "Jet"
        setattr(self, "b_n%s" % objName, tree.valueReader("n%s" % objName))
        for varName in ["eta", "phi",
                        "jetId", "puId", "btagDeepFlavB", "btagDeepFlavCvL", "btagDeepFlavCvB", "btagDeepFlavQG","hadronFlavour"]:
            setattr(self, "b_%s_%s" % (objName, varName), tree.arrayReader("%s_%s" % (objName, varName)))
        for varName in ["pt", "mass"]:
            setattr(self, "b_%s_%s" % (objName, varName), tree.arrayReader("%s_%s_%s" % (objName, varName, "nom")))

        jetVar = self.jetVar
        setattr(self, "b_Jet_Uncert", tree.arrayReader("Jet_jecUncert%s" % jetVar))

        # 1. In the __init__ method, after initializing self.worker, set the C-tag Working Point cuts.
        # The values 0.085 and 0.342 correspond to the DeepJet Medium WP for UL2016.
        self.worker.setCjetCuts(0.098, 0.370)

        self.worker.setMET(self.b_MET_pt, self.b_MET_phi)
        self.worker.setElectrons(self.b_Electron_pt, self.b_Electron_eta, self.b_Electron_isPFcand, self.b_Electron_lostHits,
                                 self.b_Electron_dxy, self.b_Electron_dz, self.b_Electron_sip3d, self.b_Electron_miniPFRelIso_all,
                                 self.b_Electron_sieie, self.b_Electron_hoe, self.b_Electron_eInvMinusPInv, self.b_Electron_convVeto,
                                 self.b_Electron_mvaTOP, self.b_Electron_phi, self.b_Electron_mass, self.b_Electron_charge,
                                 self.b_Electron_jetPt, self.b_Electron_jetEta, self.b_Electron_jetPhi, self.b_Electron_jetMass,
                                 self.b_Electron_jetPtRatio, self.b_Electron_jetRelIso, self.b_Electron_mvaFall17V2noIso_WPL,
                                 self.b_Electron_pfRelIso03_all, self.b_Electron_jetBTagDeepFlavor, self.b_Electron_deltaEtaSC)
        self.worker.setMuons(self.b_Muon_pt, self.b_Muon_eta, self.b_Muon_isPFcand, self.b_Muon_isTracker, self.b_Muon_isGlobal,
                             self.b_Muon_dxy, self.b_Muon_dz, self.b_Muon_sip3d, self.b_Muon_miniPFRelIso_all, self.b_Muon_mediumId,
                             self.b_Muon_mvaTOP, self.b_Muon_phi, self.b_Muon_mass, self.b_Muon_charge,
                             self.b_Muon_jetPt, self.b_Muon_jetEta, self.b_Muon_jetPhi, self.b_Muon_jetMass,
                             self.b_Muon_jetPtRatio, self.b_Muon_pfRelIso04_all, self.b_Muon_jetBTagDeepFlavor)
        self.worker.setJets(self.b_Jet_pt, self.b_Jet_eta, self.b_Jet_phi, self.b_Jet_mass,
                            self.b_Jet_jetId, self.b_Jet_puId, self.b_Jet_btagDeepFlavB,
                            self.b_Jet_btagDeepFlavCvL, self.b_Jet_btagDeepFlavCvB, self.b_Jet_btagDeepFlavQG,
                            self.b_Jet_hadronFlavour, self.b_Jet_Uncert)

        self._ttreereaderversion = tree._ttreereaderversion

        pass
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        if event._tree._ttreereaderversion > self._ttreereaderversion:
            self.initReaders(event._tree)
        self.worker.analyze()

        jetVar = self.jetVar

        for objName in ["ptsortGoodJet1", "ptsortGoodJet2", "ptsortGoodJet3"]:
            for varName in ["pt", "eta", "phi"]:
                setattr(event._tree, "b_out_%s_%s_%sUp" % (objName, varName, jetVar), getattr(self.worker, 'get_%s_%s' % (objName, varName))())
                self.out.fillBranch("%s_%s_%sUp" % (objName, varName, jetVar), getattr(event._tree, 'b_out_%s_%s_%sUp' % (objName, varName, jetVar)))
        for objName in ["GoodJet"]:
            for varName in ["pt", "eta", "phi", "mass"]:
                setattr(event._tree, "b_out_%s_%s_%sUp" % (objName, varName, jetVar), getattr(self.worker, 'get_%s_%s' % (objName, varName))())
                self.out.fillBranch("%s_%s_%sUp" % (objName, varName, jetVar), getattr(event._tree, 'b_out_%s_%s_%sUp' % (objName, varName, jetVar)))
        for varName in ["W_MT"]:
            setattr(event._tree, "b_out_%s_%sUp" % (varName, jetVar), getattr(self.worker, 'get_%s' % (varName))())
            self.out.fillBranch("%s_%sUp" % (varName, jetVar), getattr(event._tree, "b_out_%s_%sUp" % (varName, jetVar)))
        for varName in ["GoodJet_index", "GoodJet_DeepFlavB", "GoodJet_hadronFlavour", "nBjet"]:
            setattr(event._tree, "b_out_%s_%sUp" % (varName, jetVar), getattr(self.worker, 'get_%s' % (varName))())
            self.out.fillBranch("%s_%sUp" % (varName, jetVar), getattr(event._tree, "b_out_%s_%sUp" % (varName, jetVar)))
        ## Special care for nGoodJet
        setattr(event._tree, "b_out_nGoodJet_%sUp" % jetVar, self.worker.get_nGoodJet())
        self.out.fillBranch("nGoodJet_%sUp" % jetVar, getattr(event._tree, "b_out_nGoodJet_%sUp" % jetVar))

        correct_met_pt, correct_met_phi = self.worker.get_MET()
        setattr(event._tree, "b_out_MET_pt_%sUp" % jetVar, correct_met_pt)
        self.out.fillBranch("MET_pt_%sUp" % jetVar, getattr(event._tree, "b_out_MET_pt_%sUp" % jetVar))
        setattr(event._tree, "b_out_MET_phi_%sUp" % jetVar, correct_met_phi)
        self.out.fillBranch("MET_phi_%sUp" % jetVar, getattr(event._tree, "b_out_MET_phi_%sUp" % jetVar))

        # 3. In the analyze method, after calling the worker to process the event, fill the branches:
        nCjet = self.worker.get_nCjet()
        nUDSjet = self.worker.get_nUDSjet()
        cvl = self.worker.get_GoodJet_DeepFlavCvL()
        cvb = self.worker.get_GoodJet_DeepFlavCvB()
        qg = self.worker.get_GoodJet_DeepFlavQG()

        # Register properties in event._tree so that subsequent WeightProducer modules can read them
        setattr(event._tree, "b_out_nCjet_%sUp" % jetVar, nCjet)
        setattr(event._tree, "b_out_nUDSjet_%sUp" % jetVar, nUDSjet)
        setattr(event._tree, "b_out_GoodJet_DeepFlavCvL_%sUp" % jetVar, cvl)
        setattr(event._tree, "b_out_GoodJet_DeepFlavCvB_%sUp" % jetVar, cvb)
        setattr(event._tree, "b_out_GoodJet_DeepFlavQG_%sUp" % jetVar, qg)

        # Fill the output ROOT file with branches
        self.out.fillBranch("nCjet_%sUp" % jetVar, getattr(event._tree, "b_out_nCjet_%sUp" % jetVar))
        self.out.fillBranch("nUDSjet_%sUp" % jetVar, getattr(event._tree, "b_out_nUDSjet_%sUp" % jetVar))
        self.out.fillBranch("GoodJet_DeepFlavCvL_%sUp" % jetVar, getattr(event._tree, "b_out_GoodJet_DeepFlavCvL_%sUp" % jetVar))
        self.out.fillBranch("GoodJet_DeepFlavCvB_%sUp" % jetVar, getattr(event._tree, "b_out_GoodJet_DeepFlavCvB_%sUp" % jetVar))
        self.out.fillBranch("GoodJet_DeepFlavQG_%sUp" % jetVar, getattr(event._tree, "b_out_GoodJet_DeepFlavQG_%sUp" % jetVar))

        return True

jetUncert_fcnc_ElElEl_2016_AbsoluteStat = lambda : FCNCTriLepton(mode="ElElEl", jetVar="AbsoluteStat")
jetUncert_fcnc_MuElEl_2016_AbsoluteStat = lambda : FCNCTriLepton(mode="MuElEl", jetVar="AbsoluteStat")
jetUncert_fcnc_ElMuMu_2016_AbsoluteStat = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="AbsoluteStat")
jetUncert_fcnc_MuMuMu_2016_AbsoluteStat = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="AbsoluteStat")
jetUncert_fcnc_NPLElElEl_2016_AbsoluteStat = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="AbsoluteStat")
jetUncert_fcnc_NPLMuElEl_2016_AbsoluteStat = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="AbsoluteStat")
jetUncert_fcnc_NPLElMuMu_2016_AbsoluteStat = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="AbsoluteStat")
jetUncert_fcnc_NPLMuMuMu_2016_AbsoluteStat = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="AbsoluteStat")

jetUncert_fcnc_ElElEl_2016_AbsoluteScale = lambda : FCNCTriLepton(mode="ElElEl", jetVar="AbsoluteScale")
jetUncert_fcnc_MuElEl_2016_AbsoluteScale = lambda : FCNCTriLepton(mode="MuElEl", jetVar="AbsoluteScale")
jetUncert_fcnc_ElMuMu_2016_AbsoluteScale = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="AbsoluteScale")
jetUncert_fcnc_MuMuMu_2016_AbsoluteScale = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="AbsoluteScale")
jetUncert_fcnc_NPLElElEl_2016_AbsoluteScale = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="AbsoluteScale")
jetUncert_fcnc_NPLMuElEl_2016_AbsoluteScale = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="AbsoluteScale")
jetUncert_fcnc_NPLElMuMu_2016_AbsoluteScale = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="AbsoluteScale")
jetUncert_fcnc_NPLMuMuMu_2016_AbsoluteScale = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="AbsoluteScale")

jetUncert_fcnc_ElElEl_2016_AbsoluteMPFBias = lambda : FCNCTriLepton(mode="ElElEl", jetVar="AbsoluteMPFBias")
jetUncert_fcnc_MuElEl_2016_AbsoluteMPFBias = lambda : FCNCTriLepton(mode="MuElEl", jetVar="AbsoluteMPFBias")
jetUncert_fcnc_ElMuMu_2016_AbsoluteMPFBias = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="AbsoluteMPFBias")
jetUncert_fcnc_MuMuMu_2016_AbsoluteMPFBias = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="AbsoluteMPFBias")
jetUncert_fcnc_NPLElElEl_2016_AbsoluteMPFBias = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="AbsoluteMPFBias")
jetUncert_fcnc_NPLMuElEl_2016_AbsoluteMPFBias = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="AbsoluteMPFBias")
jetUncert_fcnc_NPLElMuMu_2016_AbsoluteMPFBias = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="AbsoluteMPFBias")
jetUncert_fcnc_NPLMuMuMu_2016_AbsoluteMPFBias = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="AbsoluteMPFBias")

jetUncert_fcnc_ElElEl_2016_Fragmentation = lambda : FCNCTriLepton(mode="ElElEl", jetVar="Fragmentation")
jetUncert_fcnc_MuElEl_2016_Fragmentation = lambda : FCNCTriLepton(mode="MuElEl", jetVar="Fragmentation")
jetUncert_fcnc_ElMuMu_2016_Fragmentation = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="Fragmentation")
jetUncert_fcnc_MuMuMu_2016_Fragmentation = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="Fragmentation")
jetUncert_fcnc_NPLElElEl_2016_Fragmentation = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="Fragmentation")
jetUncert_fcnc_NPLMuElEl_2016_Fragmentation = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="Fragmentation")
jetUncert_fcnc_NPLElMuMu_2016_Fragmentation = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="Fragmentation")
jetUncert_fcnc_NPLMuMuMu_2016_Fragmentation = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="Fragmentation")

jetUncert_fcnc_ElElEl_2016_SinglePionECAL = lambda : FCNCTriLepton(mode="ElElEl", jetVar="SinglePionECAL")
jetUncert_fcnc_MuElEl_2016_SinglePionECAL = lambda : FCNCTriLepton(mode="MuElEl", jetVar="SinglePionECAL")
jetUncert_fcnc_ElMuMu_2016_SinglePionECAL = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="SinglePionECAL")
jetUncert_fcnc_MuMuMu_2016_SinglePionECAL = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="SinglePionECAL")
jetUncert_fcnc_NPLElElEl_2016_SinglePionECAL = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="SinglePionECAL")
jetUncert_fcnc_NPLMuElEl_2016_SinglePionECAL = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="SinglePionECAL")
jetUncert_fcnc_NPLElMuMu_2016_SinglePionECAL = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="SinglePionECAL")
jetUncert_fcnc_NPLMuMuMu_2016_SinglePionECAL = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="SinglePionECAL")

jetUncert_fcnc_ElElEl_2016_SinglePionHCAL = lambda : FCNCTriLepton(mode="ElElEl", jetVar="SinglePionHCAL")
jetUncert_fcnc_MuElEl_2016_SinglePionHCAL = lambda : FCNCTriLepton(mode="MuElEl", jetVar="SinglePionHCAL")
jetUncert_fcnc_ElMuMu_2016_SinglePionHCAL = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="SinglePionHCAL")
jetUncert_fcnc_MuMuMu_2016_SinglePionHCAL = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="SinglePionHCAL")
jetUncert_fcnc_NPLElElEl_2016_SinglePionHCAL = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="SinglePionHCAL")
jetUncert_fcnc_NPLMuElEl_2016_SinglePionHCAL = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="SinglePionHCAL")
jetUncert_fcnc_NPLElMuMu_2016_SinglePionHCAL = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="SinglePionHCAL")
jetUncert_fcnc_NPLMuMuMu_2016_SinglePionHCAL = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="SinglePionHCAL")

jetUncert_fcnc_ElElEl_2016_FlavorQCD = lambda : FCNCTriLepton(mode="ElElEl", jetVar="FlavorQCD")
jetUncert_fcnc_MuElEl_2016_FlavorQCD = lambda : FCNCTriLepton(mode="MuElEl", jetVar="FlavorQCD")
jetUncert_fcnc_ElMuMu_2016_FlavorQCD = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="FlavorQCD")
jetUncert_fcnc_MuMuMu_2016_FlavorQCD = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="FlavorQCD")
jetUncert_fcnc_NPLElElEl_2016_FlavorQCD = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="FlavorQCD")
jetUncert_fcnc_NPLMuElEl_2016_FlavorQCD = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="FlavorQCD")
jetUncert_fcnc_NPLElMuMu_2016_FlavorQCD = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="FlavorQCD")
jetUncert_fcnc_NPLMuMuMu_2016_FlavorQCD = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="FlavorQCD")

jetUncert_fcnc_ElElEl_2016_FlavorPureQuark = lambda : FCNCTriLepton(mode="ElElEl", jetVar="FlavorPureQuark")
jetUncert_fcnc_MuElEl_2016_FlavorPureQuark = lambda : FCNCTriLepton(mode="MuElEl", jetVar="FlavorPureQuark")
jetUncert_fcnc_ElMuMu_2016_FlavorPureQuark = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="FlavorPureQuark")
jetUncert_fcnc_MuMuMu_2016_FlavorPureQuark = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="FlavorPureQuark")
jetUncert_fcnc_NPLElElEl_2016_FlavorPureQuark = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="FlavorPureQuark")
jetUncert_fcnc_NPLMuElEl_2016_FlavorPureQuark = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="FlavorPureQuark")
jetUncert_fcnc_NPLElMuMu_2016_FlavorPureQuark = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="FlavorPureQuark")
jetUncert_fcnc_NPLMuMuMu_2016_FlavorPureQuark = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="FlavorPureQuark")

jetUncert_fcnc_ElElEl_2016_FlavorPureCharm = lambda : FCNCTriLepton(mode="ElElEl", jetVar="FlavorPureCharm")
jetUncert_fcnc_MuElEl_2016_FlavorPureCharm = lambda : FCNCTriLepton(mode="MuElEl", jetVar="FlavorPureCharm")
jetUncert_fcnc_ElMuMu_2016_FlavorPureCharm = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="FlavorPureCharm")
jetUncert_fcnc_MuMuMu_2016_FlavorPureCharm = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="FlavorPureCharm")
jetUncert_fcnc_NPLElElEl_2016_FlavorPureCharm = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="FlavorPureCharm")
jetUncert_fcnc_NPLMuElEl_2016_FlavorPureCharm = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="FlavorPureCharm")
jetUncert_fcnc_NPLElMuMu_2016_FlavorPureCharm = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="FlavorPureCharm")
jetUncert_fcnc_NPLMuMuMu_2016_FlavorPureCharm = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="FlavorPureCharm")

jetUncert_fcnc_ElElEl_2016_FlavorPureBottom = lambda : FCNCTriLepton(mode="ElElEl", jetVar="FlavorPureBottom")
jetUncert_fcnc_MuElEl_2016_FlavorPureBottom = lambda : FCNCTriLepton(mode="MuElEl", jetVar="FlavorPureBottom")
jetUncert_fcnc_ElMuMu_2016_FlavorPureBottom = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="FlavorPureBottom")
jetUncert_fcnc_MuMuMu_2016_FlavorPureBottom = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="FlavorPureBottom")
jetUncert_fcnc_NPLElElEl_2016_FlavorPureBottom = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="FlavorPureBottom")
jetUncert_fcnc_NPLMuElEl_2016_FlavorPureBottom = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="FlavorPureBottom")
jetUncert_fcnc_NPLElMuMu_2016_FlavorPureBottom = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="FlavorPureBottom")
jetUncert_fcnc_NPLMuMuMu_2016_FlavorPureBottom = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="FlavorPureBottom")

jetUncert_fcnc_ElElEl_2016_TimePtEta = lambda : FCNCTriLepton(mode="ElElEl", jetVar="TimePtEta")
jetUncert_fcnc_MuElEl_2016_TimePtEta = lambda : FCNCTriLepton(mode="MuElEl", jetVar="TimePtEta")
jetUncert_fcnc_ElMuMu_2016_TimePtEta = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="TimePtEta")
jetUncert_fcnc_MuMuMu_2016_TimePtEta = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="TimePtEta")
jetUncert_fcnc_NPLElElEl_2016_TimePtEta = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="TimePtEta")
jetUncert_fcnc_NPLMuElEl_2016_TimePtEta = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="TimePtEta")
jetUncert_fcnc_NPLElMuMu_2016_TimePtEta = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="TimePtEta")
jetUncert_fcnc_NPLMuMuMu_2016_TimePtEta = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="TimePtEta")

jetUncert_fcnc_ElElEl_2016_RelativeJEREC1 = lambda : FCNCTriLepton(mode="ElElEl", jetVar="RelativeJEREC1")
jetUncert_fcnc_MuElEl_2016_RelativeJEREC1 = lambda : FCNCTriLepton(mode="MuElEl", jetVar="RelativeJEREC1")
jetUncert_fcnc_ElMuMu_2016_RelativeJEREC1 = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="RelativeJEREC1")
jetUncert_fcnc_MuMuMu_2016_RelativeJEREC1 = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="RelativeJEREC1")
jetUncert_fcnc_NPLElElEl_2016_RelativeJEREC1 = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="RelativeJEREC1")
jetUncert_fcnc_NPLMuElEl_2016_RelativeJEREC1 = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="RelativeJEREC1")
jetUncert_fcnc_NPLElMuMu_2016_RelativeJEREC1 = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="RelativeJEREC1")
jetUncert_fcnc_NPLMuMuMu_2016_RelativeJEREC1 = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="RelativeJEREC1")

jetUncert_fcnc_ElElEl_2016_RelativeJEREC2 = lambda : FCNCTriLepton(mode="ElElEl", jetVar="RelativeJEREC2")
jetUncert_fcnc_MuElEl_2016_RelativeJEREC2 = lambda : FCNCTriLepton(mode="MuElEl", jetVar="RelativeJEREC2")
jetUncert_fcnc_ElMuMu_2016_RelativeJEREC2 = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="RelativeJEREC2")
jetUncert_fcnc_MuMuMu_2016_RelativeJEREC2 = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="RelativeJEREC2")
jetUncert_fcnc_NPLElElEl_2016_RelativeJEREC2 = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="RelativeJEREC2")
jetUncert_fcnc_NPLMuElEl_2016_RelativeJEREC2 = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="RelativeJEREC2")
jetUncert_fcnc_NPLElMuMu_2016_RelativeJEREC2 = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="RelativeJEREC2")
jetUncert_fcnc_NPLMuMuMu_2016_RelativeJEREC2 = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="RelativeJEREC2")

jetUncert_fcnc_ElElEl_2016_RelativeJERHF = lambda : FCNCTriLepton(mode="ElElEl", jetVar="RelativeJERHF")
jetUncert_fcnc_MuElEl_2016_RelativeJERHF = lambda : FCNCTriLepton(mode="MuElEl", jetVar="RelativeJERHF")
jetUncert_fcnc_ElMuMu_2016_RelativeJERHF = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="RelativeJERHF")
jetUncert_fcnc_MuMuMu_2016_RelativeJERHF = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="RelativeJERHF")
jetUncert_fcnc_NPLElElEl_2016_RelativeJERHF = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="RelativeJERHF")
jetUncert_fcnc_NPLMuElEl_2016_RelativeJERHF = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="RelativeJERHF")
jetUncert_fcnc_NPLElMuMu_2016_RelativeJERHF = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="RelativeJERHF")
jetUncert_fcnc_NPLMuMuMu_2016_RelativeJERHF = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="RelativeJERHF")

jetUncert_fcnc_ElElEl_2016_RelativePtBB = lambda : FCNCTriLepton(mode="ElElEl", jetVar="RelativePtBB")
jetUncert_fcnc_MuElEl_2016_RelativePtBB = lambda : FCNCTriLepton(mode="MuElEl", jetVar="RelativePtBB")
jetUncert_fcnc_ElMuMu_2016_RelativePtBB = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="RelativePtBB")
jetUncert_fcnc_MuMuMu_2016_RelativePtBB = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="RelativePtBB")
jetUncert_fcnc_NPLElElEl_2016_RelativePtBB = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="RelativePtBB")
jetUncert_fcnc_NPLMuElEl_2016_RelativePtBB = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="RelativePtBB")
jetUncert_fcnc_NPLElMuMu_2016_RelativePtBB = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="RelativePtBB")
jetUncert_fcnc_NPLMuMuMu_2016_RelativePtBB = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="RelativePtBB")

jetUncert_fcnc_ElElEl_2016_RelativePtEC1 = lambda : FCNCTriLepton(mode="ElElEl", jetVar="RelativePtEC1")
jetUncert_fcnc_MuElEl_2016_RelativePtEC1 = lambda : FCNCTriLepton(mode="MuElEl", jetVar="RelativePtEC1")
jetUncert_fcnc_ElMuMu_2016_RelativePtEC1 = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="RelativePtEC1")
jetUncert_fcnc_MuMuMu_2016_RelativePtEC1 = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="RelativePtEC1")
jetUncert_fcnc_NPLElElEl_2016_RelativePtEC1 = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="RelativePtEC1")
jetUncert_fcnc_NPLMuElEl_2016_RelativePtEC1 = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="RelativePtEC1")
jetUncert_fcnc_NPLElMuMu_2016_RelativePtEC1 = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="RelativePtEC1")
jetUncert_fcnc_NPLMuMuMu_2016_RelativePtEC1 = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="RelativePtEC1")

jetUncert_fcnc_ElElEl_2016_RelativePtEC2 = lambda : FCNCTriLepton(mode="ElElEl", jetVar="RelativePtEC2")
jetUncert_fcnc_MuElEl_2016_RelativePtEC2 = lambda : FCNCTriLepton(mode="MuElEl", jetVar="RelativePtEC2")
jetUncert_fcnc_ElMuMu_2016_RelativePtEC2 = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="RelativePtEC2")
jetUncert_fcnc_MuMuMu_2016_RelativePtEC2 = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="RelativePtEC2")
jetUncert_fcnc_NPLElElEl_2016_RelativePtEC2 = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="RelativePtEC2")
jetUncert_fcnc_NPLMuElEl_2016_RelativePtEC2 = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="RelativePtEC2")
jetUncert_fcnc_NPLElMuMu_2016_RelativePtEC2 = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="RelativePtEC2")
jetUncert_fcnc_NPLMuMuMu_2016_RelativePtEC2 = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="RelativePtEC2")

jetUncert_fcnc_ElElEl_2016_RelativePtHF = lambda : FCNCTriLepton(mode="ElElEl", jetVar="RelativePtHF")
jetUncert_fcnc_MuElEl_2016_RelativePtHF = lambda : FCNCTriLepton(mode="MuElEl", jetVar="RelativePtHF")
jetUncert_fcnc_ElMuMu_2016_RelativePtHF = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="RelativePtHF")
jetUncert_fcnc_MuMuMu_2016_RelativePtHF = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="RelativePtHF")
jetUncert_fcnc_NPLElElEl_2016_RelativePtHF = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="RelativePtHF")
jetUncert_fcnc_NPLMuElEl_2016_RelativePtHF = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="RelativePtHF")
jetUncert_fcnc_NPLElMuMu_2016_RelativePtHF = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="RelativePtHF")
jetUncert_fcnc_NPLMuMuMu_2016_RelativePtHF = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="RelativePtHF")

jetUncert_fcnc_ElElEl_2016_RelativeBal = lambda : FCNCTriLepton(mode="ElElEl", jetVar="RelativeBal")
jetUncert_fcnc_MuElEl_2016_RelativeBal = lambda : FCNCTriLepton(mode="MuElEl", jetVar="RelativeBal")
jetUncert_fcnc_ElMuMu_2016_RelativeBal = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="RelativeBal")
jetUncert_fcnc_MuMuMu_2016_RelativeBal = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="RelativeBal")
jetUncert_fcnc_NPLElElEl_2016_RelativeBal = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="RelativeBal")
jetUncert_fcnc_NPLMuElEl_2016_RelativeBal = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="RelativeBal")
jetUncert_fcnc_NPLElMuMu_2016_RelativeBal = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="RelativeBal")
jetUncert_fcnc_NPLMuMuMu_2016_RelativeBal = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="RelativeBal")

jetUncert_fcnc_ElElEl_2016_RelativeSample = lambda : FCNCTriLepton(mode="ElElEl", jetVar="RelativeSample")
jetUncert_fcnc_MuElEl_2016_RelativeSample = lambda : FCNCTriLepton(mode="MuElEl", jetVar="RelativeSample")
jetUncert_fcnc_ElMuMu_2016_RelativeSample = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="RelativeSample")
jetUncert_fcnc_MuMuMu_2016_RelativeSample = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="RelativeSample")
jetUncert_fcnc_NPLElElEl_2016_RelativeSample = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="RelativeSample")
jetUncert_fcnc_NPLMuElEl_2016_RelativeSample = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="RelativeSample")
jetUncert_fcnc_NPLElMuMu_2016_RelativeSample = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="RelativeSample")
jetUncert_fcnc_NPLMuMuMu_2016_RelativeSample = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="RelativeSample")

jetUncert_fcnc_ElElEl_2016_RelativeFSR = lambda : FCNCTriLepton(mode="ElElEl", jetVar="RelativeFSR")
jetUncert_fcnc_MuElEl_2016_RelativeFSR = lambda : FCNCTriLepton(mode="MuElEl", jetVar="RelativeFSR")
jetUncert_fcnc_ElMuMu_2016_RelativeFSR = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="RelativeFSR")
jetUncert_fcnc_MuMuMu_2016_RelativeFSR = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="RelativeFSR")
jetUncert_fcnc_NPLElElEl_2016_RelativeFSR = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="RelativeFSR")
jetUncert_fcnc_NPLMuElEl_2016_RelativeFSR = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="RelativeFSR")
jetUncert_fcnc_NPLElMuMu_2016_RelativeFSR = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="RelativeFSR")
jetUncert_fcnc_NPLMuMuMu_2016_RelativeFSR = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="RelativeFSR")

jetUncert_fcnc_ElElEl_2016_RelativeStatFSR = lambda : FCNCTriLepton(mode="ElElEl", jetVar="RelativeStatFSR")
jetUncert_fcnc_MuElEl_2016_RelativeStatFSR = lambda : FCNCTriLepton(mode="MuElEl", jetVar="RelativeStatFSR")
jetUncert_fcnc_ElMuMu_2016_RelativeStatFSR = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="RelativeStatFSR")
jetUncert_fcnc_MuMuMu_2016_RelativeStatFSR = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="RelativeStatFSR")
jetUncert_fcnc_NPLElElEl_2016_RelativeStatFSR = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="RelativeStatFSR")
jetUncert_fcnc_NPLMuElEl_2016_RelativeStatFSR = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="RelativeStatFSR")
jetUncert_fcnc_NPLElMuMu_2016_RelativeStatFSR = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="RelativeStatFSR")
jetUncert_fcnc_NPLMuMuMu_2016_RelativeStatFSR = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="RelativeStatFSR")

jetUncert_fcnc_ElElEl_2016_RelativeStatEC = lambda : FCNCTriLepton(mode="ElElEl", jetVar="RelativeStatEC")
jetUncert_fcnc_MuElEl_2016_RelativeStatEC = lambda : FCNCTriLepton(mode="MuElEl", jetVar="RelativeStatEC")
jetUncert_fcnc_ElMuMu_2016_RelativeStatEC = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="RelativeStatEC")
jetUncert_fcnc_MuMuMu_2016_RelativeStatEC = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="RelativeStatEC")
jetUncert_fcnc_NPLElElEl_2016_RelativeStatEC = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="RelativeStatEC")
jetUncert_fcnc_NPLMuElEl_2016_RelativeStatEC = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="RelativeStatEC")
jetUncert_fcnc_NPLElMuMu_2016_RelativeStatEC = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="RelativeStatEC")
jetUncert_fcnc_NPLMuMuMu_2016_RelativeStatEC = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="RelativeStatEC")

jetUncert_fcnc_ElElEl_2016_RelativeStatHF = lambda : FCNCTriLepton(mode="ElElEl", jetVar="RelativeStatHF")
jetUncert_fcnc_MuElEl_2016_RelativeStatHF = lambda : FCNCTriLepton(mode="MuElEl", jetVar="RelativeStatHF")
jetUncert_fcnc_ElMuMu_2016_RelativeStatHF = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="RelativeStatHF")
jetUncert_fcnc_MuMuMu_2016_RelativeStatHF = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="RelativeStatHF")
jetUncert_fcnc_NPLElElEl_2016_RelativeStatHF = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="RelativeStatHF")
jetUncert_fcnc_NPLMuElEl_2016_RelativeStatHF = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="RelativeStatHF")
jetUncert_fcnc_NPLElMuMu_2016_RelativeStatHF = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="RelativeStatHF")
jetUncert_fcnc_NPLMuMuMu_2016_RelativeStatHF = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="RelativeStatHF")

jetUncert_fcnc_ElElEl_2016_PileUpDataMC = lambda : FCNCTriLepton(mode="ElElEl", jetVar="PileUpDataMC")
jetUncert_fcnc_MuElEl_2016_PileUpDataMC = lambda : FCNCTriLepton(mode="MuElEl", jetVar="PileUpDataMC")
jetUncert_fcnc_ElMuMu_2016_PileUpDataMC = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="PileUpDataMC")
jetUncert_fcnc_MuMuMu_2016_PileUpDataMC = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="PileUpDataMC")
jetUncert_fcnc_NPLElElEl_2016_PileUpDataMC = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="PileUpDataMC")
jetUncert_fcnc_NPLMuElEl_2016_PileUpDataMC = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="PileUpDataMC")
jetUncert_fcnc_NPLElMuMu_2016_PileUpDataMC = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="PileUpDataMC")
jetUncert_fcnc_NPLMuMuMu_2016_PileUpDataMC = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="PileUpDataMC")

jetUncert_fcnc_ElElEl_2016_PileUpPtRef = lambda : FCNCTriLepton(mode="ElElEl", jetVar="PileUpPtRef")
jetUncert_fcnc_MuElEl_2016_PileUpPtRef = lambda : FCNCTriLepton(mode="MuElEl", jetVar="PileUpPtRef")
jetUncert_fcnc_ElMuMu_2016_PileUpPtRef = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="PileUpPtRef")
jetUncert_fcnc_MuMuMu_2016_PileUpPtRef = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="PileUpPtRef")
jetUncert_fcnc_NPLElElEl_2016_PileUpPtRef = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="PileUpPtRef")
jetUncert_fcnc_NPLMuElEl_2016_PileUpPtRef = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="PileUpPtRef")
jetUncert_fcnc_NPLElMuMu_2016_PileUpPtRef = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="PileUpPtRef")
jetUncert_fcnc_NPLMuMuMu_2016_PileUpPtRef = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="PileUpPtRef")

jetUncert_fcnc_ElElEl_2016_PileUpPtBB = lambda : FCNCTriLepton(mode="ElElEl", jetVar="PileUpPtBB")
jetUncert_fcnc_MuElEl_2016_PileUpPtBB = lambda : FCNCTriLepton(mode="MuElEl", jetVar="PileUpPtBB")
jetUncert_fcnc_ElMuMu_2016_PileUpPtBB = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="PileUpPtBB")
jetUncert_fcnc_MuMuMu_2016_PileUpPtBB = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="PileUpPtBB")
jetUncert_fcnc_NPLElElEl_2016_PileUpPtBB = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="PileUpPtBB")
jetUncert_fcnc_NPLMuElEl_2016_PileUpPtBB = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="PileUpPtBB")
jetUncert_fcnc_NPLElMuMu_2016_PileUpPtBB = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="PileUpPtBB")
jetUncert_fcnc_NPLMuMuMu_2016_PileUpPtBB = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="PileUpPtBB")

jetUncert_fcnc_ElElEl_2016_PileUpPtEC1 = lambda : FCNCTriLepton(mode="ElElEl", jetVar="PileUpPtEC1")
jetUncert_fcnc_MuElEl_2016_PileUpPtEC1 = lambda : FCNCTriLepton(mode="MuElEl", jetVar="PileUpPtEC1")
jetUncert_fcnc_ElMuMu_2016_PileUpPtEC1 = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="PileUpPtEC1")
jetUncert_fcnc_MuMuMu_2016_PileUpPtEC1 = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="PileUpPtEC1")
jetUncert_fcnc_NPLElElEl_2016_PileUpPtEC1 = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="PileUpPtEC1")
jetUncert_fcnc_NPLMuElEl_2016_PileUpPtEC1 = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="PileUpPtEC1")
jetUncert_fcnc_NPLElMuMu_2016_PileUpPtEC1 = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="PileUpPtEC1")
jetUncert_fcnc_NPLMuMuMu_2016_PileUpPtEC1 = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="PileUpPtEC1")

jetUncert_fcnc_ElElEl_2016_PileUpPtEC2 = lambda : FCNCTriLepton(mode="ElElEl", jetVar="PileUpPtEC2")
jetUncert_fcnc_MuElEl_2016_PileUpPtEC2 = lambda : FCNCTriLepton(mode="MuElEl", jetVar="PileUpPtEC2")
jetUncert_fcnc_ElMuMu_2016_PileUpPtEC2 = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="PileUpPtEC2")
jetUncert_fcnc_MuMuMu_2016_PileUpPtEC2 = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="PileUpPtEC2")
jetUncert_fcnc_NPLElElEl_2016_PileUpPtEC2 = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="PileUpPtEC2")
jetUncert_fcnc_NPLMuElEl_2016_PileUpPtEC2 = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="PileUpPtEC2")
jetUncert_fcnc_NPLElMuMu_2016_PileUpPtEC2 = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="PileUpPtEC2")
jetUncert_fcnc_NPLMuMuMu_2016_PileUpPtEC2 = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="PileUpPtEC2")

jetUncert_fcnc_ElElEl_2016_PileUpPtHF = lambda : FCNCTriLepton(mode="ElElEl", jetVar="PileUpPtHF")
jetUncert_fcnc_MuElEl_2016_PileUpPtHF = lambda : FCNCTriLepton(mode="MuElEl", jetVar="PileUpPtHF")
jetUncert_fcnc_ElMuMu_2016_PileUpPtHF = lambda : FCNCTriLepton(mode="ElMuMu", jetVar="PileUpPtHF")
jetUncert_fcnc_MuMuMu_2016_PileUpPtHF = lambda : FCNCTriLepton(mode="MuMuMu", jetVar="PileUpPtHF")
jetUncert_fcnc_NPLElElEl_2016_PileUpPtHF = lambda : FCNCTriLepton(mode="NPLElElEl", jetVar="PileUpPtHF")
jetUncert_fcnc_NPLMuElEl_2016_PileUpPtHF = lambda : FCNCTriLepton(mode="NPLMuElEl", jetVar="PileUpPtHF")
jetUncert_fcnc_NPLElMuMu_2016_PileUpPtHF = lambda : FCNCTriLepton(mode="NPLElMuMu", jetVar="PileUpPtHF")
jetUncert_fcnc_NPLMuMuMu_2016_PileUpPtHF = lambda : FCNCTriLepton(mode="NPLMuMuMu", jetVar="PileUpPtHF")