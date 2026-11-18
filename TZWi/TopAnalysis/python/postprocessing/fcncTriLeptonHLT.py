import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from TZWi.TopAnalysis.postprocessing.CombineHLT import CombineHLT

###2024
setFile = "fcncTriLepton/2024.yaml"
hlt_MC2024 = lambda f=setFile: CombineHLT(fileName=f, hltSet="RunIIISummer24")
for dataset in ['EGamma', 'Muon', 'MuonEG']:
    for e in "CDEFGHI":
        vars()['hlt_Run2024%s_%s' % (e, dataset)] = lambda x=dataset, f=setFile: CombineHLT(fileName=f, hltSet="Run2024.%s" % x, doFilter=True)
