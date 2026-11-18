import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from TZWi.TopAnalysis.postprocessing.CombineHLT import CombineHLT

# Dedicated HLT factories for Run-3 2024 lepton-SF NTuple production.
# This intentionally keeps the common fcncTriLeptonHLT.py behavior untouched.
# Difference from the common module: MC is also filtered with doFilter=True so
# that data and MC enter the lepton-SF measurement in the same trigger phase space.

setFile = "fcncTriLepton/2024.yaml"

hlt_MC2024 = lambda f=setFile: CombineHLT(
    fileName=f,
    hltSet="RunIIISummer24",
    doFilter=True,
)

for dataset in ["EGamma", "Muon", "MuonEG"]:
    for era in "CDEFGHI":
        vars()["hlt_Run2024%s_%s" % (era, dataset)] = lambda x=dataset, f=setFile: CombineHLT(
            fileName=f,
            hltSet="Run2024.%s" % x,
            doFilter=True,
        )
