import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import os
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class CombineHLT(Module, object):
    def __init__(self, *args, **kwargs):
        import yaml

        self.outName = kwargs.get("outName") if "outName" in kwargs else "HLT"
        self.doFilter = kwargs.get("doFilter") if "doFilter" in kwargs else False

        if "/CombineHLTCppWorker_cc.so" not in  ROOT.gSystem.GetLibraries():
            print("Load C++ CombineHLT worker module")
            base = os.getenv("NANOAODTOOLS_BASE")
            if base:
                ROOT.gROOT.ProcessLine(".L %s/src/CombineHLTCppWorker.cc+O" % base)
            else:
                base = "%s/src/TZWi/TopAnalysis"%os.getenv("CMSSW_BASE")
                ROOT.gSystem.Load("libTZWiTopAnalysis.so")
                ROOT.gROOT.ProcessLine(".L %s/interface/CombineHLTCppWorker.h" % base)

        fName = kwargs.get("fileName")
        setName = kwargs.get("hltSet")
        base = "%s/src/TZWi/TopAnalysis"%os.getenv("CMSSW_BASE")
        d = yaml.safe_load(open(base+"/data/combineHLT/"+fName))

        # 1. Clean up the formula string and separate the operators
        formula = d[setName].replace('\n', '')
        self.oprs = ["(", ")", "||", "&&", "!"]
        for opr in self.oprs: 
            formula = formula.replace(opr, " %s " % opr)
        self.original_formula = formula.strip()

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self.outName, "O")
        self.initReaders(inputTree)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        # Delete the Worker instance first to free the C++ internal dangling pointer
        if hasattr(self, 'worker'):
            del self.worker
            
        # Clean up TTreeReader objects to prevent cppyy double-free segfaults
        for attr in list(self.__dict__.keys()):
            if attr.startswith("b_"):
                delattr(self, attr)

    def initReaders(self, tree):
        self.valid_names = []
        shortFormula = []
        
        # Safely unwrap the NanoAODTools InputTree wrapper to get the actual ROOT TTree
        actual_tree = tree._tree if hasattr(tree, '_tree') else tree

        for tok in self.original_formula.split():
            if tok in self.oprs:
                shortFormula.append(tok)
                continue

            # CRITICAL FIX: Check branch existence using the unwrapped actual TTree
            if not actual_tree.GetBranch(tok) or actual_tree.GetBranchStatus(tok) == 0:
                shortFormula.append("0")
                continue

            if tok in self.valid_names:
                idx = self.valid_names.index(tok)
            else:
                idx = len(self.valid_names)
                self.valid_names.append(tok)
            shortFormula.append("[%d]" % idx)
            
        shortFormula = ''.join(shortFormula)

        self.worker = ROOT.CombineHLTCppWorker(str(shortFormula), self.outName)
        
        for name in self.valid_names:
            setattr(self, "b_"+name, tree.valueReader(name))
            self.worker.addHLT(getattr(self, "b_"+name))

        self._ttreereaderversion = tree._ttreereaderversion

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        if event._tree._ttreereaderversion > self._ttreereaderversion:
            self.initReaders(event._tree)
            # If a new Reader is created in the middle of the loop because the file has changed,
            # the address of the current entry is remapped.
            event._tree._ttreereader.SetLocalEntry(event._tree._ttreereader.GetCurrentEntry())

        res = self.worker.analyze()
        self.out.fillBranch(self.outName, res)

        if self.doFilter: return res
        return True
