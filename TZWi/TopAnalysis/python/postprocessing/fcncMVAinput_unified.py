import ROOT
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True

import yaml
import os
from ROOT import TLorentzVector
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class FCNCMVAinput(Module, object):
    """
    Unified version of:
      - fcncMVAinput.py
      - fcncMVAinput_jetVarUp.py
      - fcncMVAinput_jetVarDown.py

    Behavior is selected only by constructor arguments / lambda names.

    Important compatibility note:
      * central, jetVarUp, and jetVarDown are all aligned to the 2018 configuration
        for xsec YAML lookup, data stream names, and integrated luminosity.
    """

    JET_VARS = [
        "AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias",
        "Fragmentation", "SinglePionECAL", "SinglePionHCAL",
        "FlavorQCD",
        "TimePtEta", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
        "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF",
        "RelativeBal", "RelativeSample", "RelativeFSR", "RelativeStatFSR",
        "RelativeStatEC", "RelativeStatHF", "PileUpDataMC", "PileUpPtRef",
        "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF",
    ]

    OUTPUT_FLOATS = [
        "MVAinput_WLZL1_dR", "MVAinput_WLZL1_dPhi", "MVAinput_WLZL1_dEta",
        "MVAinput_WLZL2_dR", "MVAinput_WLZL2_dPhi", "MVAinput_WLZL2_dEta",
        "MVAinput_ZL1ZL2_dR", "MVAinput_ZL1ZL2_dPhi", "MVAinput_ZL1ZL2_dEta",
        "MVAinput_Wlepchapt", "MVAinput_Wlepchaeta",
        "MVAinput_WJ_DeepJetB", "MVAinput_WJ_DeepJetCvL", "MVAinput_WJ_DeepJetCvB", "MVAinput_WJ_DeepJetQG", "MVAinput_WJ_pt",
        "MVAinput_WLWJ_dR", "MVAinput_WLWJ_dPhi", "MVAinput_WLWJ_dEta",
        "MVAinput_ZL1WJ_dR", "MVAinput_ZL1WJ_dPhi", "MVAinput_ZL1WJ_dEta",
        "MVAinput_ZL2WJ_dR", "MVAinput_ZL2WJ_dPhi", "MVAinput_ZL2WJ_dEta",
        "MVAinput_Z_WJ_dR", "MVAinput_Z_WJ_dPhi", "MVAinput_Z_WJ_dEta",
        "MVAinput_ZJ_DeepJetB", "MVAinput_ZJ_DeepJetCvL", "MVAinput_ZJ_DeepJetCvB", "MVAinput_ZJ_DeepJetQG", "MVAinput_ZJ_pt",
        "MVAinput_WLZJ_dR", "MVAinput_WLZJ_dPhi", "MVAinput_WLZJ_dEta",
        "MVAinput_ZL1ZJ_dR", "MVAinput_ZL1ZJ_dPhi", "MVAinput_ZL1ZJ_dEta",
        "MVAinput_ZL2ZJ_dR", "MVAinput_ZL2ZJ_dPhi", "MVAinput_ZL2ZJ_dEta",
        "MVAinput_Z_ZJ_dR", "MVAinput_Z_ZJ_dPhi", "MVAinput_Z_ZJ_dEta",
        "MVAinput_WJZJ_dR", "MVAinput_WJZJ_dPhi", "MVAinput_WJZJ_dEta",
    ]
    OUTPUT_INTS = ["MVAinput_Status"]

    def __init__(self, *args, **kwargs):
        self.jetVar = kwargs.get("jetVar")
        self.jetSyst = kwargs.get("jetSyst")
        if self.jetSyst not in (None, "Up", "Down"):
            self.jetSyst = None
        if self.jetSyst in ("Up", "Down") and not self.jetVar:
            self.jetSyst = None
            self.jetVar = None
        self._owns_xsecnorm = False

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

    def _tree_attr(self, base):
        return "b_out_%s%s" % (base, self._suffix())

    def _branch_name(self, base):
        return "%s%s" % (base, self._suffix())

    def _get_tree_value(self, tree, base, default=None):
        return getattr(tree, self._tree_attr(base), default)

    def _load_xsec_norm(self):
        info = {}
        cmssw_base = os.environ.get("CMSSW_BASE", "")
        xsecfile = "%s/src/TZWi/TopAnalysis/test/fcncTriLepton/config/crosssection.yaml" % cmssw_base

        info.update(yaml.safe_load(open(xsecfile)))

        info.update(yaml.safe_load(open(
            "%s/src/TZWi/TopAnalysis/test/fcncTriLepton/NanoAOD/2018/MC.RunIISummer18.central.yaml" % cmssw_base
        )))
        info['dataset'].update(yaml.safe_load(open(
            "%s/src/TZWi/TopAnalysis/test/fcncTriLepton/NanoAOD/2018/MC.RunIISummer18.rareprocess.yaml" % cmssw_base
        ))['dataset'])
        data_names = ["EGamma", "DoubleMuon", "SingleMuon", "MuonEG", "Tau"]
        lumi = 59800.0

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
        for items in dataitem:
            fileset0 = items[0]
            fileset1 = fileset0.split('.')[1]
            fileset.append(fileset1)
            fileroutekey0 = list(items[1].keys())
            fileroute.append(fileroutekey0)

        findsample = ""
        for j, routes in enumerate(fileroute):
            for routeitem in routes:
                if rootsample in routeitem:
                    findsample = fileset[j]

        targetxsecweight = 1.0
        targetentry = 1
        is_data = rootsample in data_names

        if (not is_data) and findsample:
            targetxsec = 1.0
            for xsecslot in xsecitem:
                if findsample == xsecslot[0]:
                    targetxsec = xsecslot[1]
            for entryslot in entryitem:
                if findsample == entryslot[0]:
                    targetentry = entryslot[1]
            targetxsecweight = (float(targetxsec) / float(targetentry)) * lumi

        self.targetxsecweight = targetxsecweight

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        # xsecNorm is a common branch shared by central / jetVarUp / jetVarDown.
        # Declare it only once, by whichever FCNCMVAinput instance is attached first.
        if not getattr(self.out, "_fcncMVAinput_xsecnorm_declared", False):
            self.out.branch("xsecNorm", "F")
            setattr(self.out, "_fcncMVAinput_xsecnorm_declared", True)
            self._owns_xsecnorm = True
        else:
            self._owns_xsecnorm = False

        for name in self.OUTPUT_FLOATS:
            self.out.branch(self._branch_name(name), "F")
        for name in self.OUTPUT_INTS:
            self.out.branch(self._branch_name(name), "i")

        self._load_xsec_norm()

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        tree = event._tree

        vals = {name: 0.0 for name in self.OUTPUT_FLOATS}
        vals["MVAinput_Status"] = 0

        if self._owns_xsecnorm:
            setattr(tree, "b_out_xsecNorm", self.targetxsecweight)
            self.out.fillBranch("xsecNorm", self.targetxsecweight)

        has_Wlep = hasattr(tree, self._tree_attr("Wlepvar"))
        has_Zlep1 = hasattr(tree, self._tree_attr("Zlep1var"))
        has_Zlep2 = hasattr(tree, self._tree_attr("Zlep2var"))
        has_leptons = has_Wlep and has_Zlep1 and has_Zlep2

        has_W_jet = hasattr(tree, self._tree_attr("jet_W"))
        has_Z_jet = hasattr(tree, self._tree_attr("jet_Z"))

        good_lepton_code = getattr(tree, "b_out_GoodLeptonCode", 0)
        nGoodJet = self._get_tree_value(tree, "nGoodJet", 0)
        nBjet = self._get_tree_value(tree, "nBjet", 0)

        if good_lepton_code == 111 and has_leptons:
            if tree.b_out_Lepton1_pdgId > 0:
                Wlepcha = -1
            else:
                Wlepcha = +1

            b_out_Wlepvar = self._get_tree_value(tree, "Wlepvar")
            b_out_Zlep1var = self._get_tree_value(tree, "Zlep1var")
            b_out_Zlep2var = self._get_tree_value(tree, "Zlep2var")

            vals["MVAinput_WLZL1_dR"] = abs(b_out_Wlepvar.DeltaR(b_out_Zlep1var))
            vals["MVAinput_WLZL1_dPhi"] = abs(b_out_Wlepvar.DeltaPhi(b_out_Zlep1var))
            vals["MVAinput_WLZL1_dEta"] = abs(b_out_Wlepvar.Eta() - b_out_Zlep1var.Eta())

            vals["MVAinput_WLZL2_dR"] = abs(b_out_Wlepvar.DeltaR(b_out_Zlep2var))
            vals["MVAinput_WLZL2_dPhi"] = abs(b_out_Wlepvar.DeltaPhi(b_out_Zlep2var))
            vals["MVAinput_WLZL2_dEta"] = abs(b_out_Wlepvar.Eta() - b_out_Zlep2var.Eta())

            vals["MVAinput_ZL1ZL2_dR"] = abs(b_out_Zlep1var.DeltaR(b_out_Zlep2var))
            vals["MVAinput_ZL1ZL2_dPhi"] = abs(b_out_Zlep1var.DeltaPhi(b_out_Zlep2var))
            vals["MVAinput_ZL1ZL2_dEta"] = abs(b_out_Zlep1var.Eta() - b_out_Zlep2var.Eta())

            vals["MVAinput_Wlepchapt"] = Wlepcha * b_out_Wlepvar.Pt()
            vals["MVAinput_Wlepchaeta"] = Wlepcha * abs(b_out_Wlepvar.Eta())

        if good_lepton_code == 111 and nGoodJet >= 1 and nBjet >= 1 and has_leptons and has_W_jet:
            b_out_Wlepvar = self._get_tree_value(tree, "Wlepvar")
            b_out_Zlep1var = self._get_tree_value(tree, "Zlep1var")
            b_out_Zlep2var = self._get_tree_value(tree, "Zlep2var")
            b_out_jet_W = self._get_tree_value(tree, "jet_W")

            vals["MVAinput_WJ_DeepJetB"] = self._get_tree_value(tree, "jet_W_DeepFlavB", 0)
            vals["MVAinput_WJ_DeepJetCvL"] = self._get_tree_value(tree, "jet_W_DeepFlavCvL", 0)
            vals["MVAinput_WJ_DeepJetCvB"] = self._get_tree_value(tree, "jet_W_DeepFlavCvB", 0)
            vals["MVAinput_WJ_DeepJetQG"] = self._get_tree_value(tree, "jet_W_DeepFlavQG", 0)
            vals["MVAinput_WJ_pt"] = b_out_jet_W.Pt()

            vals["MVAinput_WLWJ_dR"] = abs(b_out_Wlepvar.DeltaR(b_out_jet_W))
            vals["MVAinput_WLWJ_dPhi"] = abs(b_out_Wlepvar.DeltaPhi(b_out_jet_W))
            vals["MVAinput_WLWJ_dEta"] = abs(b_out_Wlepvar.Eta() - b_out_jet_W.Eta())

            vals["MVAinput_ZL1WJ_dR"] = abs(b_out_Zlep1var.DeltaR(b_out_jet_W))
            vals["MVAinput_ZL1WJ_dPhi"] = abs(b_out_Zlep1var.DeltaPhi(b_out_jet_W))
            vals["MVAinput_ZL1WJ_dEta"] = abs(b_out_Zlep1var.Eta() - b_out_jet_W.Eta())

            vals["MVAinput_ZL2WJ_dR"] = abs(b_out_Zlep2var.DeltaR(b_out_jet_W))
            vals["MVAinput_ZL2WJ_dPhi"] = abs(b_out_Zlep2var.DeltaPhi(b_out_jet_W))
            vals["MVAinput_ZL2WJ_dEta"] = abs(b_out_Zlep2var.Eta() - b_out_jet_W.Eta())

            zll = b_out_Zlep1var + b_out_Zlep2var
            vals["MVAinput_Z_WJ_dR"] = abs(zll.DeltaR(b_out_jet_W))
            vals["MVAinput_Z_WJ_dPhi"] = abs(zll.DeltaPhi(b_out_jet_W))
            vals["MVAinput_Z_WJ_dEta"] = abs(zll.Eta() - b_out_jet_W.Eta())

        if good_lepton_code == 111 and nGoodJet >= 2 and nBjet >= 1 and has_leptons and has_Z_jet:
            b_out_Wlepvar = self._get_tree_value(tree, "Wlepvar")
            b_out_Zlep1var = self._get_tree_value(tree, "Zlep1var")
            b_out_Zlep2var = self._get_tree_value(tree, "Zlep2var")
            b_out_jet_W = self._get_tree_value(tree, "jet_W")
            b_out_jet_Z = self._get_tree_value(tree, "jet_Z")

            vals["MVAinput_ZJ_DeepJetB"] = self._get_tree_value(tree, "jet_Z_DeepFlavB", 0)
            vals["MVAinput_ZJ_DeepJetCvL"] = self._get_tree_value(tree, "jet_Z_DeepFlavCvL", 0)
            vals["MVAinput_ZJ_DeepJetCvB"] = self._get_tree_value(tree, "jet_Z_DeepFlavCvB", 0)
            vals["MVAinput_ZJ_DeepJetQG"] = self._get_tree_value(tree, "jet_Z_DeepFlavQG", 0)
            vals["MVAinput_ZJ_pt"] = b_out_jet_Z.Pt()

            vals["MVAinput_WLZJ_dR"] = abs(b_out_Wlepvar.DeltaR(b_out_jet_Z))
            vals["MVAinput_WLZJ_dPhi"] = abs(b_out_Wlepvar.DeltaPhi(b_out_jet_Z))
            vals["MVAinput_WLZJ_dEta"] = abs(b_out_Wlepvar.Eta() - b_out_jet_Z.Eta())

            vals["MVAinput_ZL1ZJ_dR"] = abs(b_out_Zlep1var.DeltaR(b_out_jet_Z))
            vals["MVAinput_ZL1ZJ_dPhi"] = abs(b_out_Zlep1var.DeltaPhi(b_out_jet_Z))
            vals["MVAinput_ZL1ZJ_dEta"] = abs(b_out_Zlep1var.Eta() - b_out_jet_Z.Eta())

            vals["MVAinput_ZL2ZJ_dR"] = abs(b_out_Zlep2var.DeltaR(b_out_jet_Z))
            vals["MVAinput_ZL2ZJ_dPhi"] = abs(b_out_Zlep2var.DeltaPhi(b_out_jet_Z))
            vals["MVAinput_ZL2ZJ_dEta"] = abs(b_out_Zlep2var.Eta() - b_out_jet_Z.Eta())

            zll = b_out_Zlep1var + b_out_Zlep2var
            vals["MVAinput_Z_ZJ_dR"] = abs(zll.DeltaR(b_out_jet_Z))
            vals["MVAinput_Z_ZJ_dPhi"] = abs(zll.DeltaPhi(b_out_jet_Z))
            vals["MVAinput_Z_ZJ_dEta"] = abs(zll.Eta() - b_out_jet_Z.Eta())

            if has_W_jet:
                vals["MVAinput_WJZJ_dR"] = abs(b_out_jet_W.DeltaR(b_out_jet_Z))
                vals["MVAinput_WJZJ_dPhi"] = abs(b_out_jet_W.DeltaPhi(b_out_jet_Z))
                vals["MVAinput_WJZJ_dEta"] = abs(b_out_jet_W.Eta() - b_out_jet_Z.Eta())

        for name in self.OUTPUT_INTS:
            self.out.fillBranch(self._branch_name(name), vals[name])
        for name in self.OUTPUT_FLOATS:
            self.out.fillBranch(self._branch_name(name), vals[name])

        return True


fcncMVAinput = lambda: FCNCMVAinput()


def _make_up_lambda(jet_var):
    return lambda: FCNCMVAinput(jetVar=jet_var, jetSyst="Up")


def _make_down_lambda(jet_var):
    return lambda: FCNCMVAinput(jetVar=jet_var, jetSyst="Down")


for _jet_var in FCNCMVAinput.JET_VARS:
    globals()["fcncMVAinputUp_%s" % _jet_var] = _make_up_lambda(_jet_var)
    globals()["fcncMVAinputDown_%s" % _jet_var] = _make_down_lambda(_jet_var)
