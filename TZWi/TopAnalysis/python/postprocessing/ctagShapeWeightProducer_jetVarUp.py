import os
import correctionlib
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class CtagShapeWeightProducer(Module):
    def __init__(self, *args, **kwargs):
        self.jetVar = kwargs.get("jetVar")
        cmssw_base = os.environ.get("CMSSW_BASE", "")
        # Updated to use the new ctagging.json.gz file
        self.json_path = os.path.join(cmssw_base, "src/TZWi/TopAnalysis/data/BTV/UL2018/ctagging.json.gz")
        
        if not os.path.exists(self.json_path):
            raise RuntimeError("C-tagging JSON file not found: %s" % self.json_path)
            
        self.systs = [
            "Extrap", "Interp", "Stat", 
            "XSec_BRUnc_DYJets_b", "XSec_BRUnc_DYJets_c", "XSec_BRUnc_WJets_c",
            "jer", "jesTotal", "puWeight"
        ]
        
        print(f"Loading C-tagging shape SF from correctionlib: {self.json_path}")
        
        # Load correctionlib
        self.cset = correctionlib.CorrectionSet.from_file(self.json_path)
        # BTV POG standard key for c-tagging reshaping is typically "deepJet_shape" or "deepJet_cshape"
        # We assume "deepJet_shape" as it is standard, but you can change it if your JSON uses a different key.
        self.ctag_corr = self.cset["deepJet_shape"]

    def evaluate_corr(self, corr, **kwargs):
        """
        Dynamically maps keyword arguments to the exact positional arguments 
        required by the specific correctionlib JSON schema.
        Includes a fallback to 'central' (or 1.0) for missing systematics,
        emulating the old CSV parser's safe fallback behavior.
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
                args.append(kwargs.get("working_point", "shape"))
            else:
                args.append(kwargs.get(inp.name, 0.0))
                
        try:
            # Try to evaluate the exact systematic and flavor combination
            return corr.evaluate(*args)
        except Exception:
            # Catch correctionlib IndexError if the systematic is not available for this flavor.
            # Fallback to the 'central' scale factor, just like the old CSV parser did.
            if syst_index >= 0 and args[syst_index] != "central":
                args[syst_index] = "central"
                try:
                    return corr.evaluate(*args)
                except Exception:
                    # If even 'central' fails (e.g., completely undefined bin), return 1.0
                    return 1.0
            return 1.0

    def beginJob(self): pass

    def endJob(self):
        # Explicitly delete correctionlib objects to prevent pybind11/ROOT GC segfaults at exit
        if hasattr(self, 'ctag_corr'): del self.ctag_corr
        if hasattr(self, 'cset'): del self.cset

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        jetVar = self.jetVar
        
        self.out.branch("CtagWeight_%sUp" % jetVar, "F")
        for syst in self.systs:
            self.out.branch("CtagWeight_up_%s_%sUp" % (syst, jetVar), "F")
            self.out.branch("CtagWeight_down_%s_%sUp" % (syst, jetVar), "F")

    def analyze(self, event):
        jetVar = self.jetVar
        weight_central = 1.0
        weights_up = {syst: 1.0 for syst in self.systs}
        weights_down = {syst: 1.0 for syst in self.systs}

        # Safely bypass the event getattr issue
        n_good_jet = getattr(event._tree, "b_out_nGoodJet_%sUp" % jetVar, 0)
        
        if n_good_jet > 0:
            pts = getattr(event._tree, "b_out_GoodJet_pt_%sUp" % jetVar)
            etas = getattr(event._tree, "b_out_GoodJet_eta_%sUp" % jetVar)
            flavs = getattr(event._tree, "b_out_GoodJet_hadronFlavour_%sUp" % jetVar)
            cvls = getattr(event._tree, "b_out_GoodJet_DeepFlavCvL_%sUp" % jetVar)
            cvbs = getattr(event._tree, "b_out_GoodJet_DeepFlavCvB_%sUp" % jetVar)

            for jet_index in range(n_good_jet):
                jet_pt = pts[jet_index]
                jet_eta = etas[jet_index]
                
                if jet_pt < 25.0 or abs(jet_eta) > 2.4:
                    continue

                flav = abs(flavs[jet_index])
                flav_map = 0
                if flav == 5:
                    flav_map = 5
                elif flav == 4:
                    flav_map = 4

                cvl = float(cvls[jet_index])
                cvb = float(cvbs[jet_index])

                # Prepare the full dictionary of possible variables
                eval_kwargs = {
                    "systematic": "central",
                    "flavor": flav_map,
                    "cvl": cvl,
                    "cvb": cvb,
                    "pt": jet_pt,
                    "eta": abs(jet_eta),
                    "abseta": abs(jet_eta)
                }

                sf_central = self.evaluate_corr(self.ctag_corr, **eval_kwargs)
                weight_central *= sf_central

                for syst in self.systs:
                    eval_kwargs["systematic"] = "up_" + syst
                    sf_up = self.evaluate_corr(self.ctag_corr, **eval_kwargs)
                    
                    eval_kwargs["systematic"] = "down_" + syst
                    sf_down = self.evaluate_corr(self.ctag_corr, **eval_kwargs)
                        
                    weights_up[syst] *= sf_up
                    weights_down[syst] *= sf_down

        self.out.fillBranch("CtagWeight_%sUp" % jetVar, weight_central)
        for syst in self.systs:
            self.out.fillBranch("CtagWeight_up_%s_%sUp" % (syst, jetVar), weights_up[syst])
            self.out.fillBranch("CtagWeight_down_%s_%sUp" % (syst, jetVar), weights_down[syst])

        return True


ctagWeight_AbsoluteStat = lambda : CtagShapeWeightProducer(jetVar="AbsoluteStat")
ctagWeight_AbsoluteScale = lambda : CtagShapeWeightProducer(jetVar="AbsoluteScale")
ctagWeight_AbsoluteMPFBias = lambda : CtagShapeWeightProducer(jetVar="AbsoluteMPFBias")
ctagWeight_Fragmentation = lambda : CtagShapeWeightProducer(jetVar="Fragmentation")
ctagWeight_SinglePionECAL = lambda : CtagShapeWeightProducer(jetVar="SinglePionECAL")
ctagWeight_SinglePionHCAL = lambda : CtagShapeWeightProducer(jetVar="SinglePionHCAL")
ctagWeight_FlavorQCD = lambda : CtagShapeWeightProducer(jetVar="FlavorQCD")
ctagWeight_FlavorPureQuark = lambda : CtagShapeWeightProducer(jetVar="FlavorPureQuark")
ctagWeight_FlavorPureCharm = lambda : CtagShapeWeightProducer(jetVar="FlavorPureCharm")
ctagWeight_FlavorPureBottom = lambda : CtagShapeWeightProducer(jetVar="FlavorPureBottom")
ctagWeight_TimePtEta = lambda : CtagShapeWeightProducer(jetVar="TimePtEta")
ctagWeight_RelativeJEREC1 = lambda : CtagShapeWeightProducer(jetVar="RelativeJEREC1")
ctagWeight_RelativeJEREC2 = lambda : CtagShapeWeightProducer(jetVar="RelativeJEREC2")
ctagWeight_RelativeJERHF = lambda : CtagShapeWeightProducer(jetVar="RelativeJERHF")
ctagWeight_RelativePtBB = lambda : CtagShapeWeightProducer(jetVar="RelativePtBB")
ctagWeight_RelativePtEC1 = lambda : CtagShapeWeightProducer(jetVar="RelativePtEC1")
ctagWeight_RelativePtEC2 = lambda : CtagShapeWeightProducer(jetVar="RelativePtEC2")
ctagWeight_RelativePtHF = lambda : CtagShapeWeightProducer(jetVar="RelativePtHF")
ctagWeight_RelativeBal = lambda : CtagShapeWeightProducer(jetVar="RelativeBal")
ctagWeight_RelativeSample = lambda : CtagShapeWeightProducer(jetVar="RelativeSample")
ctagWeight_RelativeFSR = lambda : CtagShapeWeightProducer(jetVar="RelativeFSR")
ctagWeight_RelativeStatFSR = lambda : CtagShapeWeightProducer(jetVar="RelativeStatFSR")
ctagWeight_RelativeStatEC = lambda : CtagShapeWeightProducer(jetVar="RelativeStatEC")
ctagWeight_RelativeStatHF = lambda : CtagShapeWeightProducer(jetVar="RelativeStatHF")
ctagWeight_PileUpDataMC = lambda : CtagShapeWeightProducer(jetVar="PileUpDataMC")
ctagWeight_PileUpPtRef = lambda : CtagShapeWeightProducer(jetVar="PileUpPtRef")
ctagWeight_PileUpPtBB = lambda : CtagShapeWeightProducer(jetVar="PileUpPtBB")
ctagWeight_PileUpPtEC1 = lambda : CtagShapeWeightProducer(jetVar="PileUpPtEC1")
ctagWeight_PileUpPtEC2 = lambda : CtagShapeWeightProducer(jetVar="PileUpPtEC2")
ctagWeight_PileUpPtHF = lambda : CtagShapeWeightProducer(jetVar="PileUpPtHF")
