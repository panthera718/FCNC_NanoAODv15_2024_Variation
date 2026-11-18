#!/bin/bash
set -euo pipefail

OUTDIR=${OUTDIR:-fakerate_2024}
INPUT_BASE=${INPUT_BASE:-ntuple_2024_FR/reco}
mkdir -p "$OUTDIR"

python3 make_fake_rate_promptMVA_2024.py \
  --flavor el \
  --input "$INPUT_BASE/FREl/JetMET*/*.root" \
  --out "$OUTDIR/fakerate_electron_promptMVA_2024_JetMET.root" \
  --csv "$OUTDIR/fakerate_electron_promptMVA_2024_JetMET.csv"

python3 make_fake_rate_promptMVA_2024.py \
  --flavor mu \
  --input "$INPUT_BASE/FRMu/JetMET*/*.root" \
  --out "$OUTDIR/fakerate_muon_promptMVA_2024_JetMET.root" \
  --csv "$OUTDIR/fakerate_muon_promptMVA_2024_JetMET.csv"
