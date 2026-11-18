import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from TZWi.TopAnalysis.postprocessing.CombineHLT import CombineHLT

# Put the YAML file in:
#   $CMSSW_BASE/src/TZWi/TopAnalysis/data/combineHLT/fcncTriLepton/2024_triggerSF.yaml
setFile = "fcncTriLepton/2024_triggerSF.yaml"

# MC modules.  These do not filter events; they only write boolean branches.
trigSFAnalysis_MC2024 = lambda f=setFile: CombineHLT(
    outName="HLT_analysis", fileName=f, hltSet="RunIIISummer24.Analysis", doFilter=False
)
trigSFReference_MC2024 = lambda f=setFile: CombineHLT(
    outName="HLT_reference", fileName=f, hltSet="RunIIISummer24.Reference", doFilter=False
)

# Optional diagnostic groups.
trigSFElectron_MC2024 = lambda f=setFile: CombineHLT(outName="HLT_group_e", fileName=f, hltSet="TriggerGroup.Electron", doFilter=False)
trigSFDoubleElectron_MC2024 = lambda f=setFile: CombineHLT(outName="HLT_group_ee", fileName=f, hltSet="TriggerGroup.DoubleElectron", doFilter=False)
trigSFTripleElectron_MC2024 = lambda f=setFile: CombineHLT(outName="HLT_group_eee", fileName=f, hltSet="TriggerGroup.TripleElectron", doFilter=False)
trigSFMuon_MC2024 = lambda f=setFile: CombineHLT(outName="HLT_group_mu", fileName=f, hltSet="TriggerGroup.Muon", doFilter=False)
trigSFDoubleMuon_MC2024 = lambda f=setFile: CombineHLT(outName="HLT_group_mumu", fileName=f, hltSet="TriggerGroup.DoubleMuon", doFilter=False)
trigSFTripleMuon_MC2024 = lambda f=setFile: CombineHLT(outName="HLT_group_mumumu", fileName=f, hltSet="TriggerGroup.TripleMuon", doFilter=False)
trigSFEMu_MC2024 = lambda f=setFile: CombineHLT(outName="HLT_group_emu", fileName=f, hltSet="TriggerGroup.EMu", doFilter=False)
trigSFMixed3L_MC2024 = lambda f=setFile: CombineHLT(outName="HLT_group_mixed3l", fileName=f, hltSet="TriggerGroup.Mixed3L", doFilter=False)

# Data modules.  HLT_DATASET in the shell script is used only to find this
# factory name.  The formulas themselves are deliberately non-exclusive.  PD
# overlap removal is not part of a trigger-efficiency measurement with an
# independent reference dataset.
for dataset in [
    "EGamma", "Muon", "MuonEG",
    "MET", "JetMET", "MET0", "MET1", "JetMET0", "JetMET1",
]:
    for era in "CDEFGHI":
        vars()["trigSFAnalysis_Run2024%s_%s" % (era, dataset)] = (
            lambda f=setFile: CombineHLT(
                outName="HLT_analysis", fileName=f, hltSet="Run2024.Analysis", doFilter=False
            )
        )
        vars()["trigSFReference_Run2024%s_%s" % (era, dataset)] = (
            lambda f=setFile: CombineHLT(
                outName="HLT_reference", fileName=f, hltSet="Run2024.Reference", doFilter=False
            )
        )
        vars()["trigSFElectron_Run2024%s_%s" % (era, dataset)] = (
            lambda f=setFile: CombineHLT(outName="HLT_group_e", fileName=f, hltSet="TriggerGroup.Electron", doFilter=False)
        )
        vars()["trigSFDoubleElectron_Run2024%s_%s" % (era, dataset)] = (
            lambda f=setFile: CombineHLT(outName="HLT_group_ee", fileName=f, hltSet="TriggerGroup.DoubleElectron", doFilter=False)
        )
        vars()["trigSFTripleElectron_Run2024%s_%s" % (era, dataset)] = (
            lambda f=setFile: CombineHLT(outName="HLT_group_eee", fileName=f, hltSet="TriggerGroup.TripleElectron", doFilter=False)
        )
        vars()["trigSFMuon_Run2024%s_%s" % (era, dataset)] = (
            lambda f=setFile: CombineHLT(outName="HLT_group_mu", fileName=f, hltSet="TriggerGroup.Muon", doFilter=False)
        )
        vars()["trigSFDoubleMuon_Run2024%s_%s" % (era, dataset)] = (
            lambda f=setFile: CombineHLT(outName="HLT_group_mumu", fileName=f, hltSet="TriggerGroup.DoubleMuon", doFilter=False)
        )
        vars()["trigSFTripleMuon_Run2024%s_%s" % (era, dataset)] = (
            lambda f=setFile: CombineHLT(outName="HLT_group_mumumu", fileName=f, hltSet="TriggerGroup.TripleMuon", doFilter=False)
        )
        vars()["trigSFEMu_Run2024%s_%s" % (era, dataset)] = (
            lambda f=setFile: CombineHLT(outName="HLT_group_emu", fileName=f, hltSet="TriggerGroup.EMu", doFilter=False)
        )
        vars()["trigSFMixed3L_Run2024%s_%s" % (era, dataset)] = (
            lambda f=setFile: CombineHLT(outName="HLT_group_mixed3l", fileName=f, hltSet="TriggerGroup.Mixed3L", doFilter=False)
        )
