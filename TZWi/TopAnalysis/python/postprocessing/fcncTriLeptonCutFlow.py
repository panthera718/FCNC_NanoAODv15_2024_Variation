import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import os
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class FCNCTriLeptonCutFlow(Module, object):
    def __init__(self, *args, **kwargs):
        #super(FCNCTriLepton, self).__init__(*args, **kwargs)
        self.mode = kwargs.get("mode")
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("CutStep", "i")

        pass
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        cutStep = 0
        while True:

            if abs(event._tree.b_out_GoodLeptonCode) != 111: break
            if ( self.mode[:3]!="NPL" ):
                if ( event._tree.b_out_nGoodLepton != 3 ): break
            if ( self.mode[:3]=="NPL" ):
                if ( event._tree.b_out_nNPLepton != 3 ): break
                fake1 = fake2 = fake3 = 0
                if event._tree.b_out_Lepton1_isFake: fake1 = 1
                if event._tree.b_out_Lepton2_isFake: fake2 = 1
                if event._tree.b_out_Lepton3_isFake: fake3 = 1
                if ( fake1 + fake2 + fake3 == 0 ): break

            if ( event._tree.b_out_ptsortLepton1_pt < 25 or event._tree.b_out_ptsortLepton2_pt < 20 or event._tree.b_out_ptsortLepton3_pt < 20 ): break
            if ( event._tree.b_out_Z_charge != 0 ): break
            if ( abs(event._tree.b_out_Z_mass-91.2) > 15 ): break
#            if ( event._tree.b_out_JetVetoSum > 9 ): break
#            if ( event._tree.b_out_nGoodJet < 1 ): break
#            if ( event._tree.b_out_nBjet < 1 ): break

            cutStep += 1

            break

        self.out.fillBranch("CutStep", cutStep)

        return cutStep > 0

cutFlow_MuMuMu = lambda: FCNCTriLeptonCutFlow(mode="MuMuMu")
cutFlow_ElElEl = lambda: FCNCTriLeptonCutFlow(mode="ElElEl")
cutFlow_MuElEl = lambda: FCNCTriLeptonCutFlow(mode="MuElEl")
cutFlow_ElMuMu = lambda: FCNCTriLeptonCutFlow(mode="ElMuMu")
cutFlow_NPLMuMuMu = lambda: FCNCTriLeptonCutFlow(mode="NPLMuMuMu")
cutFlow_NPLElElEl = lambda: FCNCTriLeptonCutFlow(mode="NPLElElEl")
cutFlow_NPLMuElEl = lambda: FCNCTriLeptonCutFlow(mode="NPLMuElEl")
cutFlow_NPLElMuMu = lambda: FCNCTriLeptonCutFlow(mode="NPLElMuMu")

cutFlow_TaTaTa = lambda: FCNCTriLeptonCutFlow(mode="TaTaTa")
cutFlow_ElTaTa = lambda: FCNCTriLeptonCutFlow(mode="ElTaTa")
cutFlow_MuTaTa = lambda: FCNCTriLeptonCutFlow(mode="MuTaTa")
cutFlow_TaElEl = lambda: FCNCTriLeptonCutFlow(mode="TaElEl")
cutFlow_TaMuMu = lambda: FCNCTriLeptonCutFlow(mode="TaMuMu")
cutFlow_NPLTaTaTa = lambda: FCNCTriLeptonCutFlow(mode="NPLTaTaTa")
cutFlow_NPLElTaTa = lambda: FCNCTriLeptonCutFlow(mode="NPLElTaTa")
cutFlow_NPLMuTaTa = lambda: FCNCTriLeptonCutFlow(mode="NPLMuTaTa")
cutFlow_NPLTaElEl = lambda: FCNCTriLeptonCutFlow(mode="NPLTaElEl")
cutFlow_NPLTaMuMu = lambda: FCNCTriLeptonCutFlow(mode="NPLTaMuMu")
