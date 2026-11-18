#!/bin/bash
set -euo pipefail

OUTDIR=${OUTDIR:-fakerate_2024}
INPUT_BASE=${INPUT_BASE:-ntuple_2024_FR/reco}
ETA_MODE=${ETA_MODE:-signed}
PYTHON=${PYTHON:-python3}
FR_MAKER=${FR_MAKER:-make_fake_rate_promptMVA_2024_signedEta.py}
BTAG_BRANCH=${BTAG_BRANCH:-FR_nBJet}

# Nominal recommendation:
#   inclusive: reference map
#   0b       : FR_nBJet == 0
#   ge1b     : FR_nBJet >= 1
#
# To additionally produce 1b and >=2b maps:
#   MAKE_DETAILED_BTAG=1 ./run_make_fake_rates_2024_signedEta.sh
#
# Or explicitly choose categories:
#   BTAG_CATEGORIES="inclusive 0b 1b ge2b" ./run_make_fake_rates_2024_signedEta.sh
BTAG_CATEGORIES=${BTAG_CATEGORIES:-"inclusive 0b ge1b"}
if [[ "${MAKE_DETAILED_BTAG:-0}" == "1" ]]; then
  BTAG_CATEGORIES="$BTAG_CATEGORIES 1b ge2b"
fi

mkdir -p "$OUTDIR"

case "$ETA_MODE" in
  signed) ETA_SUFFIX="SignedEta" ;;
  abs)    ETA_SUFFIX="AbsEta" ;;
  *)
    echo "[run_make_fake_rates_2024_signedEta.sh] ERROR: ETA_MODE must be signed or abs, got '$ETA_MODE'" >&2
    exit 2
    ;;
esac

category_suffix() {
  case "$1" in
    inclusive) echo "Inclusive" ;;
    0b)        echo "0B" ;;
    ge1b)      echo "GE1B" ;;
    1b)        echo "1B" ;;
    ge2b)      echo "GE2B" ;;
    *)
      echo "[run_make_fake_rates_2024_signedEta.sh] ERROR: unsupported b-tag category '$1'" >&2
      return 2
      ;;
  esac
}

run_one() {
  local flavor=$1
  local flavor_name=$2
  local mode_dir=$3
  local category=$4
  local btag_suffix
  btag_suffix=$(category_suffix "$category")

  local stem="fakerate_${flavor_name}_promptMVA_2024_JetMET_${ETA_SUFFIX}_${btag_suffix}"

  echo
  echo "[run_make_fake_rates_2024_signedEta.sh] flavor=$flavor category=$category"
  echo "[run_make_fake_rates_2024_signedEta.sh] output=$OUTDIR/${stem}.root"

  "$PYTHON" "$FR_MAKER" \
    --flavor "$flavor" \
    --eta-mode "$ETA_MODE" \
    --btag-category "$category" \
    --btag-branch "$BTAG_BRANCH" \
    --input "$INPUT_BASE/$mode_dir/JetMET*/*.root" \
    --out "$OUTDIR/${stem}.root" \
    --csv "$OUTDIR/${stem}.csv"
}

# Avoid duplicate work when MAKE_DETAILED_BTAG adds an already requested category.
declare -A SEEN=()
for category in $BTAG_CATEGORIES; do
  if [[ -n "${SEEN[$category]:-}" ]]; then
    continue
  fi
  SEEN[$category]=1

  run_one el electron FREl "$category"
  run_one mu muon FRMu "$category"
done

echo
echo "[run_make_fake_rates_2024_signedEta.sh] completed categories: ${!SEEN[*]}"
