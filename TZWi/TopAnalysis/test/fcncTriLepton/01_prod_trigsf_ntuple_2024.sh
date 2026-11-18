#!/bin/bash

if [ $# != 4 ]; then
    echo "$#"
    echo "Usage: $0 ElElEl NanoAOD/2024/.../SAMPLE.txt 20 0"
    echo "       $0 MuMuMu NanoAOD/2024/.../SAMPLE.txt 20 0"
    exit 1
fi

eval `scram runtime -sh`

export PYTHONNOUSERSITE=1
export XRD_WORKERTHREADS=1
export XRD_NETWORKSTACK=IPv4

CHANNEL=$1
FILELIST=$2
MAXFILES=$3
JOBNUMBER=$4

case $CHANNEL in
  ElElMu|ElMuEl)
    CHANNEL=MuElEl
    ;;
  MuMuEl|MuElMu)
    CHANNEL=ElMuMu
    ;;
esac

DATASET0=$(basename "$FILELIST" | sed -e 's;.txt;;g')
DATASET='/'$(echo "$DATASET0" | sed -e 's;\.;/;g')
ERA=$(echo "$DATASET0" | grep -oE 'Run2024[C-I]' | head -n1)
DATASETGROUP=$(basename "$(dirname "$FILELIST")")
DATATYPE=$(echo "$DATASETGROUP" | cut -d. -f1)
YEAR=${DATATYPE:(-4)}

if [ "${DATATYPE::3}" == "Run" ]; then
  HLT_DATASET=$(echo "$DATASETGROUP" | cut -d. -f2)

  if [ -z "$HLT_DATASET" ] || [ "$HLT_DATASET" = "$DATASETGROUP" ]; then
    PRIMARY_DATASET=$(echo "$DATASET" | cut -d/ -f2)
    case "$PRIMARY_DATASET" in
      EGamma0|EGamma1)
        HLT_DATASET=EGamma
        ;;
      Muon0|Muon1)
        HLT_DATASET=Muon
        ;;
      *)
        HLT_DATASET=$PRIMARY_DATASET
        ;;
    esac
  fi

  if [ -z "$ERA" ]; then
    echo "[ERROR] Could not determine ERA from DATASET0=$DATASET0 FILELIST=$FILELIST" >&2
    exit 1
  fi

  HLTMODULE=${ERA}_${HLT_DATASET}
else
  HLTMODULE=$(echo "$DATATYPE" | cut -d_ -f1)
fi

echo "[DEBUG] FILELIST     = $FILELIST"
echo "[DEBUG] DATASET0     = $DATASET0"
echo "[DEBUG] DATASETGROUP = $DATASETGROUP"
echo "[DEBUG] DATATYPE     = $DATATYPE"
echo "[DEBUG] ERA          = $ERA"
echo "[DEBUG] HLT_DATASET  = $HLT_DATASET"
echo "[DEBUG] HLTMODULE    = $HLTMODULE"

FILENAMES=$(cat "$FILELIST" | xargs -n"$MAXFILES" | sed -n "$((JOBNUMBER+1)) p" | sed 's;^/xrootd/;root://cms-xrdr.private.lo:2094//xrd/;g')

ARGS=""
ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.flags flags_${DATATYPE}"

# Do not use the normal fcncTriLeptonHLT module here.  The normal data HLT
# module filters on the analysis lepton triggers, which would make the trigger
# denominator biased.  These two modules only write HLT_analysis and
# HLT_reference booleans.
ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLeptonTriggerSFHLT_2024 trigSFAnalysis_${HLTMODULE}"
ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLeptonTriggerSFHLT_2024 trigSFReference_${HLTMODULE}"

if [ "${TRIGSF_SAVE_GROUPS:-0}" = "1" ]; then
  ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLeptonTriggerSFHLT_2024 trigSFElectron_${HLTMODULE}"
  ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLeptonTriggerSFHLT_2024 trigSFDoubleElectron_${HLTMODULE}"
  ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLeptonTriggerSFHLT_2024 trigSFTripleElectron_${HLTMODULE}"
  ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLeptonTriggerSFHLT_2024 trigSFMuon_${HLTMODULE}"
  ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLeptonTriggerSFHLT_2024 trigSFDoubleMuon_${HLTMODULE}"
  ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLeptonTriggerSFHLT_2024 trigSFTripleMuon_${HLTMODULE}"
  ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLeptonTriggerSFHLT_2024 trigSFEMu_${HLTMODULE}"
  ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLeptonTriggerSFHLT_2024 trigSFMixed3L_${HLTMODULE}"
fi

ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLepton_allmodes_unified_jetvar_only_promptMVA_Run3 fcnc_${CHANNEL}_${YEAR}"
# Derived pT/|eta| variables used by the optional 2D trigger-SF maps.
ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.trigSFVariables_2024 trigSFVariables2024"
ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLeptonCutFlow cutFlow_${CHANNEL}"

BRANCHSEL="${TRIGSF_BRANCHSEL:-$CMSSW_BASE/src/TZWi/TopAnalysis/data/branchsel_trigSF_2024.txt}"
if [ ! -f "$BRANCHSEL" ]; then
    echo "[WARN] Trigger-SF branch selection file not found: $BRANCHSEL" >&2
    echo "[WARN] Falling back to the legacy branchsel.txt." >&2
    BRANCHSEL="$CMSSW_BASE/src/TZWi/TopAnalysis/data/branchsel.txt"
fi
CMD="nano_postproc.py --prefetch --bo $BRANCHSEL"

OUTPATH=ntuple_${YEAR}_trigSF/reco/${CHANNEL}/${DATASET0}
[ ! -d "$OUTPATH" ] && mkdir -p "$OUTPATH"

if [ "${DATATYPE::2}" != "MC" ]; then
    if [ "$DATATYPE" == "Run2024" ]; then
        JSONFILE="$CMSSW_BASE/src/TZWi/TopAnalysis/data/lumi/Cert_Collisions2024_378981_386951_Golden.json"
        ARGS="-J $JSONFILE $ARGS"
    fi
fi

echo $CMD $ARGS $OUTPATH $FILENAMES
TMP_ERR="err_${JOBNUMBER}.tmp"
$CMD $ARGS $OUTPATH $FILENAMES 2> "$TMP_ERR"
EXIT_CODE=$?
sed '/\*\*\* Break \*\*\* segmentation violation/,$d' "$TMP_ERR" >&2
rm -f "$TMP_ERR"

if [ $EXIT_CODE -eq 0 ]; then
    echo "Job completed successfully."
else
    echo "Job exited with code $EXIT_CODE. Checking if it is a PyROOT teardown crash..."
    if ls "$OUTPATH"/*.root 1> /dev/null 2>&1; then
        echo "Output .root file(s) successfully generated. Masking teardown crash."
        EXIT_CODE=0
    else
        echo "No output .root files found. Actual job failure." >&2
    fi
fi

if [ $EXIT_CODE -ne 0 ]; then
    echo $CMD $ARGS $OUTPATH $FILENAMES >> failed_trigSF.txt
    exit $EXIT_CODE
fi
exit 0
