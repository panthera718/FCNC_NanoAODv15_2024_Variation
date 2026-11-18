#!/usr/bin/env python3
import glob
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

patterns = sys.argv[1:] or [
    "ntuple_2024_FR/reco/FREl/JetMET*/*.root",
    "ntuple_2024_FR/reco/FRMu/JetMET*/*.root",
]
files = []
for p in patterns:
    files.extend(sorted(glob.glob(p)))
files = list(dict.fromkeys(files))
print("[check_fr_ntuple] nfiles =", len(files))

branches = [
    "FRRegion", "FRTrigger_pass", "FRProbe_flavor", "FRProbe_pt",
    "FRProbe_eta", "FRProbe_absEta", "FRProbe_isLoose", "FRProbe_isTight",
    "FRProbe_leadJetPt", "FRProbe_met_plus_mt",
]

total = 0
pos_eta = 0
neg_eta = 0
zero_eta = 0
for path in files[:20]:
    tf = ROOT.TFile.Open(path)
    tree = tf.Get("Events") if tf else None
    entries = tree.GetEntries() if tree else -1
    total += max(entries, 0)
    print("\n", path)
    print("  entries =", entries)
    if tree:
        for br in branches:
            print("  %-24s %s" % (br, "OK" if tree.GetBranch(br) else "MISSING"))
        if tree.GetBranch("FRProbe_eta"):
            ncheck = min(entries, 20000)
            for i in range(ncheck):
                tree.GetEntry(i)
                eta = float(tree.FRProbe_eta)
                if eta > 0.0:
                    pos_eta += 1
                elif eta < 0.0:
                    neg_eta += 1
                else:
                    zero_eta += 1
    if tf:
        tf.Close()
print("\n[check_fr_ntuple] displayed files =", min(len(files), 20))
print("[check_fr_ntuple] entries in displayed files =", total)
print("[check_fr_ntuple] sampled eta sign counts: eta>0=%d eta<0=%d eta=0=%d" % (pos_eta, neg_eta, zero_eta))
