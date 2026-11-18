import os
import correctionlib
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

def is_relevant_syst_for_shape_corr(flavor_btv, syst):
    """
    Returns true if a flavor/syst combination is relevant.
    Note: correctionlib uses standard PDG ID flavors: 5(b), 4(c), 0(light)
    """
    if flavor_btv == 5:   # b-jet
        return syst in ["jes", "lf", "hfstats1", "hfstats2"]
    elif flavor_btv == 4: # c-jet
        return syst in ["cferr1", "cferr2"]
    elif flavor_btv == 0: # light-jet
        return syst in ["jes", "hf", "lfstats1", "lfstats2"]
    return False

class BtagShapeWeightProducer(Module):
    def __init__(self):
        cmssw_base = os.environ.get("CMSSW_BASE", "")
        
        # Update path to use the json.gz file instead of the CSV
        self.json_path = os.path.join(cmssw_base, "src/TZWi/TopAnalysis/data/BTV/UL2018/btagging.json.gz")
        
        if not os.path.exists(self.json_path):
            raise RuntimeError("BTV JSON file not found: %s" % self.json_path)

        print(f"Loading BTagCalibration from correctionlib: {self.json_path}")
        
        # Load the CorrectionSet
        self.cset = correctionlib.CorrectionSet.from_file(self.json_path)
        
        # The key for reshaping in the BTV POG standard JSON is typically "deepJet_shape"
        self.btag_corr = self.cset["deepJet_shape"]
        
        self.systs = [
            "jes", "lf", "hf", 
            "hfstats1", "hfstats2", 
            "lfstats1", "lfstats2", 
            "cferr1", "cferr2"
        ]

    def beginJob(self): pass

    def endJob(self):
        # Explicitly delete correctionlib objects to prevent pybind11/ROOT GC segfaults at exit
        if hasattr(self, 'btag_corr'): del self.btag_corr
        if hasattr(self, 'cset'): del self.cset

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        self.out.branch("BtagWeight", "F")
        for syst in self.systs:
            self.out.branch("BtagWeight_up_%s" % syst, "F")
            self.out.branch("BtagWeight_down_%s" % syst, "F")

    def analyze(self, event):
        weight_central = 1.0
        weights_up = {syst: 1.0 for syst in self.systs}
        weights_down = {syst: 1.0 for syst in self.systs}

        # FIX: Read dynamically generated variables directly from event._tree
        n_good_jet = getattr(event._tree, "b_out_nGoodJet", 0)
        
        if n_good_jet > 0:
            # FIX: Also read the arrays from event._tree
            pts = event._tree.b_out_GoodJet_pt
            etas = event._tree.b_out_GoodJet_eta
            flavs = event._tree.b_out_GoodJet_hadronFlavour
            scores = event._tree.b_out_GoodJet_DeepFlavB

            for jet_index in range(n_good_jet):
                jet_pt = pts[jet_index]
                jet_eta = etas[jet_index]
                
                if jet_pt < 25.0 or abs(jet_eta) > 2.4:
                    continue

                # Map flavors according to the JSON format: 5(b), 4(c), 0(light)
                flav = abs(flavs[jet_index])
                flav_map = 0 # default to light
                if flav == 5:
                    flav_map = 5
                elif flav == 4:
                    flav_map = 4

                score = float(scores[jet_index])

                # Evaluate using correctionlib
                sf_central = self.btag_corr.evaluate("central", flav_map, abs(jet_eta), jet_pt, score)

                if sf_central < 0.01: 
                    sf_central = 1.0
                    
                weight_central *= sf_central

                for syst in self.systs:
                    if is_relevant_syst_for_shape_corr(flav_map, syst):
                        sf_up = self.btag_corr.evaluate("up_" + syst, flav_map, abs(jet_eta), jet_pt, score)
                        sf_down = self.btag_corr.evaluate("down_" + syst, flav_map, abs(jet_eta), jet_pt, score)

                        if sf_up < 0.01: sf_up = 1.0
                        if sf_down < 0.01: sf_down = 1.0
                    else:
                        # Apply Central SF for unrelated Systematics
                        sf_up = sf_central
                        sf_down = sf_central
                        
                    weights_up[syst] *= sf_up
                    weights_down[syst] *= sf_down

        self.out.fillBranch("BtagWeight", weight_central)
        for syst in self.systs:
            self.out.fillBranch("BtagWeight_up_%s" % syst, weights_up[syst])
            self.out.fillBranch("BtagWeight_down_%s" % syst, weights_down[syst])

        return True

btagWeight = lambda: BtagShapeWeightProducer()
