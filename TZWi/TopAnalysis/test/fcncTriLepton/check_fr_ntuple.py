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
for pattern in patterns:
    files.extend(sorted(glob.glob(pattern)))
files = list(dict.fromkeys(files))
if not files:
    raise SystemExit("[check_fr_ntuple_btag] no files matched")

required = [
    "FRRegion",
    "FRProbe_flavor",
    "FRProbe_pt",
    "FRProbe_eta",
    "FR_nAnalysisJet",
    "FR_nBJet",
    "FR_btagCategory0GE1",
    "FR_btagCategory01GE2",
    "FR_btagWPCode",
    "FR_btagWPValue",
    "FR_maxBtagUParTAK4B",
    "FR_leadBJetPt",
    "FRProbe_assocJetUParTAK4B",
    "FRProbe_assocJetIsBTagged",
]

chain = ROOT.TChain("Events")
for path in files:
    chain.Add(path)

missing = [name for name in required if not chain.GetBranch(name)]
if missing:
    raise SystemExit("[check_fr_ntuple_btag] missing branches: " + ", ".join(missing))

counts = {"0b": 0, "1b": 0, "ge2b": 0}
bad_category = 0
wp_values = set()
for event in chain:
    n_bjet = int(event.FR_nBJet)
    cat2 = int(event.FR_btagCategory0GE1)
    cat3 = int(event.FR_btagCategory01GE2)
    wp_values.add(round(float(event.FR_btagWPValue), 6))

    if n_bjet == 0:
        counts["0b"] += 1
    elif n_bjet == 1:
        counts["1b"] += 1
    else:
        counts["ge2b"] += 1

    if cat2 != (0 if n_bjet == 0 else 1):
        bad_category += 1
    if cat3 != min(n_bjet, 2):
        bad_category += 1

print("[check_fr_ntuple_btag] files =", len(files))
print("[check_fr_ntuple_btag] entries =", chain.GetEntries())
print("[check_fr_ntuple_btag] nBJet categories =", counts)
print("[check_fr_ntuple_btag] WP values =", sorted(wp_values))
print("[check_fr_ntuple_btag] category inconsistencies =", bad_category)

if bad_category:
    raise SystemExit(2)
