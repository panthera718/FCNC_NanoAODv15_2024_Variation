from TZWi.TopAnalysis.postprocessing.CombineHLT import CombineHLT

flags_MC2024 = lambda : CombineHLT(outName="Flag", fileName="flags/2024.yaml", hltSet="RunIIISummer24", doFilter=True)
flags_Run2024 = lambda : CombineHLT(outName="Flag", fileName="flags/2024.yaml", hltSet="Run2024", doFilter=True)
