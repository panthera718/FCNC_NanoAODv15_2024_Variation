import os
import correctionlib
import ROOT
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class QGtagWeightProducer(Module):
    def __init__(self, *args, **kwargs):
        self.jetVar = kwargs.get("jetVar")
        cmssw_base = os.environ.get("CMSSW_BASE", "")
        self.json_path = os.path.join(cmssw_base, "src/TZWi/TopAnalysis/data/BTV/UL2018/qgtagging_v0.json")
        
        if not os.path.exists(self.json_path):
            raise RuntimeError("Q/G-tagging JSON file not found: %s" % self.json_path)
            
        print(f"Loading Q/G-tagging SF from correctionlib: {self.json_path}")
        
        self.cset = correctionlib.CorrectionSet.from_file(self.json_path)
        # Q/G tagging is usually keyed as "deepJet_shape" or "deepJet_qg" depending on the JSON
        # Update this key if your JSON uses a different name for QG tagging
        self.qg_corr = self.cset["deepJet_shape"]
                        
    def evaluate_corr(self, corr, **kwargs):
        """
        Dynamically maps keyword arguments to the exact positional arguments.
        Includes a fallback to 'central' (or 1.0) for missing systematics.
        """
        args = []
        syst_index = -1
        
        for i, inp in enumerate(corr.inputs):
            name = inp.name.lower()
            if name == "systematic":
                syst_index = i
                
            if name in kwargs:
                args.append(kwargs[name])
            elif name == "working_point":
                args.append(kwargs.get("working_point", "M"))
            elif name == "flavor" or name == "absflavor":
                # Ensure flavor is always an int if missed in kwargs
                args.append(int(kwargs.get(name, 0)))
            else:
                args.append(kwargs.get(inp.name, 0.0))
                
        try:
            return corr.evaluate(*args)
        except Exception:
            # Fallback to central if a specific systematic is missing for this bin
            if syst_index >= 0 and args[syst_index] != "central":
                args[syst_index] = "central"
                try:
                    return corr.evaluate(*args)
                except Exception:
                    return 1.0
            return 1.0

    def beginJob(self): pass

    def endJob(self):
        # Explicitly delete correctionlib objects to prevent pybind11/ROOT GC segfaults at exit
        if hasattr(self, 'qg_corr'): del self.qg_corr
        if hasattr(self, 'cset'): del self.cset

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        jetVar = self.jetVar
        
        self.out.branch("QGtagWeight_%sDown" % jetVar, "F")
        self.out.branch("QGtagWeight_up_%sDown" % jetVar, "F")
        self.out.branch("QGtagWeight_down_%sDown" % jetVar, "F")

    def analyze(self, event):
        jetVar = self.jetVar
        weight_central = 1.0
        weight_up = 1.0
        weight_down = 1.0

        n_good_jet = getattr(event._tree, "b_out_nGoodJet_%sDown" % jetVar, 0)

        if n_good_jet > 0:
            pts = getattr(event._tree, "b_out_GoodJet_pt_%sDown" % jetVar)
            etas = getattr(event._tree, "b_out_GoodJet_eta_%sDown" % jetVar)
            orig_idxs = getattr(event._tree, "b_out_GoodJet_index_%sDown" % jetVar)
            qg_scores = getattr(event._tree, "b_out_GoodJet_DeepFlavQG_%sDown" % jetVar)

            for jet_index in range(n_good_jet):
                jet_pt = pts[jet_index]
                jet_eta = etas[jet_index]
                abseta = abs(jet_eta)
                
                # QG tagging SF is generally applied to the interval pt > 30, |eta| < 2.5
                if jet_pt < 30.0 or abseta > 2.5:
                    continue

                # Safely retrieve parton flavor, force casting to int
                orig_jet_idx = orig_idxs[jet_index]
                try:
                    parton_flavor = int(abs(event.Jet_partonFlavour[orig_jet_idx]))
                except (AttributeError, RuntimeError, IndexError):
                    parton_flavor = int(abs(getattr(event._tree, "b_out_GoodJet_hadronFlavour_%sDown" % jetVar)[jet_index]))

                # "Flipped likelihood" calculation: 1.0 - DeepJet
                qvg = float(1.0 - float(qg_scores[jet_index]))

                # Prepare kwargs. Explicitly provide both 'flavor' and 'absflavor' as INT types
                eval_kwargs = {
                    "systematic": "central",
                    "pt": float(jet_pt),
                    "eta": float(abseta),
                    "abseta": float(abseta),
                    "flavor": parton_flavor,
                    "absflavor": parton_flavor,
                    "qvg": qvg,
                    "discr": qvg
                }

                sf_central = self.evaluate_corr(self.qg_corr, **eval_kwargs)
                
                eval_kwargs["systematic"] = "up"
                sf_up = self.evaluate_corr(self.qg_corr, **eval_kwargs)
                
                eval_kwargs["systematic"] = "down"
                sf_down = self.evaluate_corr(self.qg_corr, **eval_kwargs)
                
                weight_central *= sf_central
                weight_up *= sf_up
                weight_down *= sf_down

        self.out.fillBranch("QGtagWeight_%sDown" % jetVar, weight_central)
        self.out.fillBranch("QGtagWeight_up_%sDown" % jetVar, weight_up)
        self.out.fillBranch("QGtagWeight_down_%sDown" % jetVar, weight_down)

        return True


qgtagWeight_AbsoluteStat = lambda : QGtagWeightProducer(jetVar="AbsoluteStat")
qgtagWeight_AbsoluteScale = lambda : QGtagWeightProducer(jetVar="AbsoluteScale")
qgtagWeight_AbsoluteMPFBias = lambda : QGtagWeightProducer(jetVar="AbsoluteMPFBias")
qgtagWeight_Fragmentation = lambda : QGtagWeightProducer(jetVar="Fragmentation")
qgtagWeight_SinglePionECAL = lambda : QGtagWeightProducer(jetVar="SinglePionECAL")
qgtagWeight_SinglePionHCAL = lambda : QGtagWeightProducer(jetVar="SinglePionHCAL")
qgtagWeight_FlavorQCD = lambda : QGtagWeightProducer(jetVar="FlavorQCD")
qgtagWeight_FlavorPureQuark = lambda : QGtagWeightProducer(jetVar="FlavorPureQuark")
qgtagWeight_FlavorPureCharm = lambda : QGtagWeightProducer(jetVar="FlavorPureCharm")
qgtagWeight_FlavorPureBottom = lambda : QGtagWeightProducer(jetVar="FlavorPureBottom")
qgtagWeight_TimePtEta = lambda : QGtagWeightProducer(jetVar="TimePtEta")
qgtagWeight_RelativeJEREC1 = lambda : QGtagWeightProducer(jetVar="RelativeJEREC1")
qgtagWeight_RelativeJEREC2 = lambda : QGtagWeightProducer(jetVar="RelativeJEREC2")
qgtagWeight_RelativeJERHF = lambda : QGtagWeightProducer(jetVar="RelativeJERHF")
qgtagWeight_RelativePtBB = lambda : QGtagWeightProducer(jetVar="RelativePtBB")
qgtagWeight_RelativePtEC1 = lambda : QGtagWeightProducer(jetVar="RelativePtEC1")
qgtagWeight_RelativePtEC2 = lambda : QGtagWeightProducer(jetVar="RelativePtEC2")
qgtagWeight_RelativePtHF = lambda : QGtagWeightProducer(jetVar="RelativePtHF")
qgtagWeight_RelativeBal = lambda : QGtagWeightProducer(jetVar="RelativeBal")
qgtagWeight_RelativeSample = lambda : QGtagWeightProducer(jetVar="RelativeSample")
qgtagWeight_RelativeFSR = lambda : QGtagWeightProducer(jetVar="RelativeFSR")
qgtagWeight_RelativeStatFSR = lambda : QGtagWeightProducer(jetVar="RelativeStatFSR")
qgtagWeight_RelativeStatEC = lambda : QGtagWeightProducer(jetVar="RelativeStatEC")
qgtagWeight_RelativeStatHF = lambda : QGtagWeightProducer(jetVar="RelativeStatHF")
qgtagWeight_PileUpDataMC = lambda : QGtagWeightProducer(jetVar="PileUpDataMC")
qgtagWeight_PileUpPtRef = lambda : QGtagWeightProducer(jetVar="PileUpPtRef")
qgtagWeight_PileUpPtBB = lambda : QGtagWeightProducer(jetVar="PileUpPtBB")
qgtagWeight_PileUpPtEC1 = lambda : QGtagWeightProducer(jetVar="PileUpPtEC1")
qgtagWeight_PileUpPtEC2 = lambda : QGtagWeightProducer(jetVar="PileUpPtEC2")
qgtagWeight_PileUpPtHF = lambda : QGtagWeightProducer(jetVar="PileUpPtHF")
