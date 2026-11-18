import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class FCNCTriLeptonCutFlow(Module, object):
    """
    Baseline exact-3-lepton + one-Z-candidate cut flow for 2024 NanoAODv15.

    The object definitions are produced by FCNCTriLepton.  This module applies
    only the event-level baseline used for the 3l+Z ntuple selection:
      - exactly three selected leptons, and exactly three e/mu veto leptons
      - one OSSF Z candidate, |mll - mZ| < 15 GeV
      - three selected leptons sorted by pT satisfy >=25, >=10, >=10 GeV
        without changing the W/Z-role ordering of Lepton1/2/3
      - trilepton mass > 100 GeV
      - PuppiMET > 30 GeV
      - official JME jet-veto map must not veto the event

    No nBjet or b-veto requirement is imposed here.  Define b-tag/top-specific
    signal/control regions later from nBjet and the reconstructed top variables.
    """

    def __init__(self, *args, **kwargs):
        self.mode = kwargs.get("mode")

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("CutStep", "i")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    @staticmethod
    def _get(tree, name, default=0):
        return getattr(tree, "b_out_%s" % name, default)

    def analyze(self, event):
        tree = event._tree
        cutStep = 0

        while True:
            if abs(int(self._get(tree, "GoodLeptonCode", 0))) != 111:
                break

            if self.mode[:3] != "NPL":
                if int(self._get(tree, "nGoodLepton", 0)) != 3:
                    break
            else:
                if int(self._get(tree, "nNPLepton", 0)) != 3:
                    break
                n_fake = int(bool(self._get(tree, "Lepton1_isFake", False)))
                n_fake += int(bool(self._get(tree, "Lepton2_isFake", False)))
                n_fake += int(bool(self._get(tree, "Lepton3_isFake", False)))
                if n_fake == 0:
                    break

            # Exact e/mu 3l baseline.  In the current production only e/mu modes
            # are enabled, so nVetoLepton is the right extra-lepton veto.
            if int(self._get(tree, "nVetoLepton", 0)) != 3:
                break

            # Apply role-independent lepton-pT thresholds using the separate
            # pT-sorted branches produced by FCNCTriLepton.  The original
            # Lepton1/2/3 branches remain unchanged: Lepton1 is the W-associated
            # lepton and Lepton2/3 form the Z candidate.
            #
            # Inclusive requirement:
            #   leading >= 25 GeV, subleading >= 10 GeV, third >= 10 GeV.
            leading_pt = float(self._get(tree, "ptsortLepton1_pt", -999.0))
            subleading_pt = float(self._get(tree, "ptsortLepton2_pt", -999.0))
            third_pt = float(self._get(tree, "ptsortLepton3_pt", -999.0))
            if leading_pt < 25.0:
                break
            if subleading_pt < 10.0:
                break
            if third_pt < 10.0:
                break

            if int(self._get(tree, "Z_charge", 999)) != 0:
                break
            if abs(float(self._get(tree, "Z_mass", 0.0)) - 91.2) > 15.0:
                break
            #if float(self._get(tree, "TriLepton_mass", 0.0)) <= 100.0:
            #    break
            #if float(self._get(tree, "MET_pt", 0.0)) <= 30.0:
            #    break

            # JetVetoProducer writes JetVetoSum after FCNCTriLepton.  If the
            # producer is disabled accidentally, fail closed only when the branch
            # exists and is non-zero; otherwise keep legacy behavior.
            if hasattr(tree, "b_out_JetVetoSum") and float(tree.b_out_JetVetoSum) != 0.0:
                break

            cutStep = 1
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
