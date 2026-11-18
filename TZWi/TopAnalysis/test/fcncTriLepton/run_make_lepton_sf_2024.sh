#!/bin/bash
set -euo pipefail

# Edit this to the SF NTuple base directory if needed.
export LEPSF_BASE=${LEPSF_BASE:-ntuple_2024_leptonSF}
CONFIG=${CONFIG:-leptonSF_samples_2024.yaml}
PYTHON=${PYTHON:-python3}
SCRIPT=${SCRIPT:-make_lepton_sf_promptMVA.py}

export LEPSF_MODE=LepSFEl
${PYTHON} ${SCRIPT} \
  --config "${CONFIG}" \
  --flavor el \
  --out leptonSF_electron_2024_signedEta.root

export LEPSF_MODE=LepSFMu
${PYTHON} ${SCRIPT} \
  --config "${CONFIG}" \
  --flavor mu \
  --out leptonSF_muon_2024_signedEta.root
