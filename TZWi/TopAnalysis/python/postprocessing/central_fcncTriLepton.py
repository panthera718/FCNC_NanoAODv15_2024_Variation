import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import os
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class FCNCTriLepton(Module, object):
    def __init__(self, *args, **kwargs):
        #super(FCNCTriLepton, self).__init__(*args, **kwargs)
        self.mode = kwargs.get("mode")
#        self.eleIdName = kwargs.get("eleId") if "eleId" in kwargs else "cutBased"
#        self.eleVidMap = kwargs.get("eleVidMap") if "eleVidMap" in kwargs else "vidNestedWPBitmap"
#        self.doNonPromptLepton = kwargs.get("doNonPromptLepton") if "doNonPromptLepton" in kwargs else False

        if "/ForBTag_FCNCTriLeptonCppWorker_cc.so" not in  ROOT.gSystem.GetLibraries():
            print("Load C++ FCNCTriLepton worker module")
            base = os.getenv("NANOAODTOOLS_BASE")
            if base:
                ROOT.gROOT.ProcessLine(".L %s/src/ForBTag_FCNCTriLeptonCppWorker.cc+O" % base)
            else:
                base = "%s/src/TZWi/TopAnalysis"%os.getenv("CMSSW_BASE")
#                ROOT.gSystem.Load("libPhysicsToolsNanoAODTools.so")
                ROOT.gSystem.Load("libTZWiTopAnalysis.so")
                ROOT.gROOT.ProcessLine(".L %s/interface/ForBTag_FCNCTriLeptonCppWorker.h" % base)

        pass
    def beginJob(self):
#        self.worker = ROOT.ForBTag_FCNCTriLeptonCppWorker(self.mode, self.doNonPromptLepton)
        self.worker = ROOT.ForBTag_FCNCTriLeptonCppWorker(self.mode)
        pass
#    def endJob(self):
#        pass

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
        for objName in ["Lepton1", "Lepton2", "Lepton3"]:
            for varName in ["pt", "eta", "phi", "mass", "PFIso", "jetPtRatio", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "jetBTagDeepFlavor"]:
                self.out.branch("%s_%s" % (objName, varName), "F")
        for objName in ["ptsortLepton1", "ptsortLepton2", "ptsortLepton3"]:
            for varName in ["pt", "eta", "phi"]:
                self.out.branch("%s_%s" % (objName, varName), "F")
        for objName in ["Z", "LeadingElectron", "LeadingMuon"]:
            for varName in ["pt", "eta", "phi", "mass"]:
                self.out.branch("%s_%s" % (objName, varName), "F")
        self.out.branch("MET_pt", "F")
        self.out.branch("MET_phi", "F")
        self.out.branch("Lepton1_pdgId", "I")
        self.out.branch("Lepton2_pdgId", "I")
        self.out.branch("Lepton3_pdgId", "I")
        self.out.branch("Lepton1_isTight", "O")
        self.out.branch("Lepton2_isTight", "O")
        self.out.branch("Lepton3_isTight", "O")
        self.out.branch("Lepton1_isFake", "O")
        self.out.branch("Lepton2_isFake", "O")
        self.out.branch("Lepton3_isFake", "O")
        self.out.branch("Lepton1_mvaTOP", "F")
        self.out.branch("Lepton2_mvaTOP", "F")
        self.out.branch("Lepton3_mvaTOP", "F")
        self.out.branch("LeadingLepton_pt", "F")
        self.out.branch("TriLepton_mass", "F")
        self.out.branch("TriLepton_pt", "F")
        self.out.branch("TriLepton_WleptonZdPhi", "F")
        self.out.branch("TriLepton_WleptonZdR", "F")
        self.out.branch("Z_charge", "I")
        self.out.branch("nVetoLepton", "i")
        self.out.branch("nGoodLepton", "i")
        self.out.branch("nNPLepton", "i")
        self.out.branch("GoodLeptonCode", "I")
        self.out.branch("nGoodElectron", "i")
        self.out.branch("nGoodMuon", "i")
        #self.out.branch("nGoodJet", "i")
        self.out.branch("GoodJet_index", "i", lenVar="nGoodJet")
        for varName in ["pt", "eta", "phi", "mass", "DeepFlavB"]:
            self.out.branch("GoodJet_%s" % varName, "F", lenVar="nGoodJet")
        for varName in ["hadronFlavour"]:
            self.out.branch("GoodJet_%s" % varName, "I", lenVar="nGoodJet")
        for objName in ["ptsortGoodJet1", "ptsortGoodJet2", "ptsortGoodJet3"]:
            for varName in ["pt", "eta", "phi"]:
                self.out.branch("%s_%s" % (objName, varName), "F")
        self.out.branch("nBjet", "i")
        self.out.branch("W_MT", "F")

        # 2. In the beginFile method, define the new branches:
        self.out.branch("nCjet", "i")
        self.out.branch("nUDSjet", "i")
        self.out.branch("GoodJet_DeepFlavCvL", "F", lenVar="nGoodJet")
        self.out.branch("GoodJet_DeepFlavCvB", "F", lenVar="nGoodJet")
        self.out.branch("GoodJet_DeepFlavQG", "F", lenVar="nGoodJet")

        self.initReaders(inputTree)
        pass
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def initReaders(self,tree):

#        self.b_MET_pt = tree.valueReader("MET_pt")
#        self.b_MET_phi = tree.valueReader("MET_phi")

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

#        objName = "Jet"
#        setattr(self, "b_n%s" % objName, tree.valueReader("n%s" % objName))
#        for varName in ["pt", "eta", "phi", "mass",
#                        "jetId", "puId", "btagDeepFlavB",]:
#            setattr(self, "b_%s_%s" % (objName, varName), tree.arrayReader("%s_%s" % (objName, varName)))

        objName = "Jet"
        setattr(self, "b_n%s" % objName, tree.valueReader("n%s" % objName))
        for varName in ["eta", "phi",
                        "jetId", "puId", "btagDeepFlavB", "btagDeepFlavCvL", "btagDeepFlavCvB", "btagDeepFlavQG","hadronFlavour"]:
            setattr(self, "b_%s_%s" % (objName, varName), tree.arrayReader("%s_%s" % (objName, varName)))
        for varName in ["pt", "mass"]:
            setattr(self, "b_%s_%s" % (objName, varName), tree.arrayReader("%s_%s_%s" % (objName, varName, "nom")))

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
                            self.b_Jet_hadronFlavour)

        self._ttreereaderversion = tree._ttreereaderversion

        pass
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        if event._tree._ttreereaderversion > self._ttreereaderversion:
            self.initReaders(event._tree)
        self.worker.analyze()

        for objName in ["Lepton1", "Lepton2", "Lepton3"]:
            for varName in ["pt", "eta", "phi", "mass", "PFIso", "jetPtRatio",
                            "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "jetBTagDeepFlavor",
                            "pdgId", "isTight", "isFake", "mvaTOP"]:
                setattr(event._tree, "b_out_%s_%s" % (objName, varName), getattr(self.worker, 'get_%s_%s' % (objName, varName))())
                self.out.fillBranch("%s_%s" % (objName, varName), getattr(event._tree, 'b_out_%s_%s' % (objName, varName)))
        for objName in ["ptsortLepton1", "ptsortLepton2", "ptsortLepton3"]:
            for varName in ["pt", "eta", "phi"]:
                setattr(event._tree, "b_out_%s_%s" % (objName, varName), getattr(self.worker, 'get_%s_%s' % (objName, varName))())
                self.out.fillBranch("%s_%s" % (objName, varName), getattr(event._tree, 'b_out_%s_%s' % (objName, varName)))
        for objName in ["ptsortGoodJet1", "ptsortGoodJet2", "ptsortGoodJet3"]:
            for varName in ["pt", "eta", "phi"]:
                setattr(event._tree, "b_out_%s_%s" % (objName, varName), getattr(self.worker, 'get_%s_%s' % (objName, varName))())
                self.out.fillBranch("%s_%s" % (objName, varName), getattr(event._tree, 'b_out_%s_%s' % (objName, varName)))
        for objName in ["Z", "GoodJet", "LeadingElectron", "LeadingMuon"]:
            for varName in ["pt", "eta", "phi", "mass"]:
                setattr(event._tree, "b_out_%s_%s" % (objName, varName), getattr(self.worker, 'get_%s_%s' % (objName, varName))())
                self.out.fillBranch("%s_%s" % (objName, varName), getattr(event._tree, 'b_out_%s_%s' % (objName, varName)))
        for varName in ["MET_pt", "MET_phi", "LeadingLepton_pt",
                        "TriLepton_mass", "TriLepton_pt", "TriLepton_WleptonZdPhi", "TriLepton_WleptonZdR",
                        "nVetoLepton", "nGoodElectron", "nGoodMuon", "nGoodLepton", "nNPLepton", "GoodLeptonCode", "Z_charge", "W_MT",
                        #"nGoodJet", #We do not keep nGoodJet here, it have to be done by the framework
                        "GoodJet_index", "GoodJet_DeepFlavB", "GoodJet_hadronFlavour",
                        "nBjet"]:
            setattr(event._tree, "b_out_%s" % (varName), getattr(self.worker, 'get_%s' % (varName))())
            self.out.fillBranch(varName, getattr(event._tree, "b_out_%s" % varName))
        ## Special care for nGoodJet
        setattr(event._tree, "b_out_nGoodJet", self.worker.get_nGoodJet())

        # 3. In the analyze method, after calling the worker to process the event, fill the branches:
        nCjet = self.worker.get_nCjet()
        nUDSjet = self.worker.get_nUDSjet()
        cvl = self.worker.get_GoodJet_DeepFlavCvL()
        cvb = self.worker.get_GoodJet_DeepFlavCvB()
        qg = self.worker.get_GoodJet_DeepFlavQG()

        # Register properties in event._tree so that subsequent WeightProducer modules can read them
        setattr(event._tree, "b_out_nCjet", nCjet)
        setattr(event._tree, "b_out_nUDSjet", nUDSjet)
        setattr(event._tree, "b_out_GoodJet_DeepFlavCvL", cvl)
        setattr(event._tree, "b_out_GoodJet_DeepFlavCvB", cvb)
        setattr(event._tree, "b_out_GoodJet_DeepFlavQG", qg)

        # Fill the output ROOT file with branches
        self.out.fillBranch("nCjet", nCjet)
        self.out.fillBranch("nUDSjet", nUDSjet)
        self.out.fillBranch("GoodJet_DeepFlavCvL", cvl)
        self.out.fillBranch("GoodJet_DeepFlavCvB", cvb)
        self.out.fillBranch("GoodJet_DeepFlavQG", qg)

        return True

central_fcnc_MuMuMu_2016 = lambda : FCNCTriLepton(mode="MuMuMu")
central_fcnc_ElElEl_2016 = lambda : FCNCTriLepton(mode="ElElEl")
central_fcnc_ElMuMu_2016 = lambda : FCNCTriLepton(mode="ElMuMu")
central_fcnc_MuElEl_2016 = lambda : FCNCTriLepton(mode="MuElEl")
central_fcnc_NPLMuMuMu_2016 = lambda : FCNCTriLepton(mode="NPLMuMuMu")
central_fcnc_NPLElElEl_2016 = lambda : FCNCTriLepton(mode="NPLElElEl")
central_fcnc_NPLElMuMu_2016 = lambda : FCNCTriLepton(mode="NPLElMuMu")
central_fcnc_NPLMuElEl_2016 = lambda : FCNCTriLepton(mode="NPLMuElEl")
