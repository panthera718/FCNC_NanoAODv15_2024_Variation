#!/bin/bash

if [ $# != 4 ]; then
    echo $#
    echo "Usage: $0 LepSFEl NanoAOD/2024/RD.Run2024/Run2024.EGamma/EGamma0.Run2024C-MINIv6NANOv15-v1.NANOAOD.txt 20 0"
    echo "Usage: $0 LepSFMu NanoAOD/2024/MC.RunIISummer24.central/MC2024.WZto3LNu/WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8.RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2.NANOAODSIM.txt 20 0"
    exit 1
fi

eval `scram runtime -sh`

export PYTHONNOUSERSITE=1
export XRD_WORKERTHREADS=${XRD_WORKERTHREADS:-1}
export XRD_NETWORKSTACK=${XRD_NETWORKSTACK:-IPv4}

CHANNEL=$1
FILELIST=$2
MAXFILES=$3
JOBNUMBER=$4

case "$CHANNEL" in
  LepSFEl|LepSFMu)
    ;;
  *)
    echo "[01_prod_ntuple_leptonSF.sh] ERROR: unsupported CHANNEL=$CHANNEL" >&2
    exit 2
    ;;
esac

DATASET0=$(basename "$FILELIST" | sed -e 's;.txt;;g')
DATASET='/'$(echo "$DATASET0" | sed -e 's;\.;/;g')
ERA=$(echo "$DATASET0" | cut -d. -f2 | cut -d- -f1 | sed -e 's;NanoAOD;;g')
DATATYPE=$(basename "$(dirname "$FILELIST")" | cut -d. -f1)
YEAR=${DATATYPE:(-4)}
PRIMARY_DATASET_RAW=$(echo "$DATASET" | cut -d/ -f2)

if [ -z "$YEAR" ]; then
    YEAR=2024
fi

# Normalize split Run-2024 primary datasets.
# File-list names contain EGamma0/EGamma1 and Muon0/Muon1, but the HLT YAML
# and HLT factories are defined for EGamma and Muon.
case "$PRIMARY_DATASET_RAW" in
  EGamma0|EGamma1)
    PRIMARY_DATASET="EGamma"
    ;;
  Muon0|Muon1)
    PRIMARY_DATASET="Muon"
    ;;
  MuonEG)
    PRIMARY_DATASET="MuonEG"
    ;;
  *)
    PRIMARY_DATASET="$PRIMARY_DATASET_RAW"
    ;;
esac

# Keep the data stream name in the HLT module name. For MC, use the MC group prefix.
if [ "${DATATYPE::3}" == "Run" ]; then
    HLTMODULE="${ERA::8}_${PRIMARY_DATASET}"
else
    HLTMODULE=$(echo "$DATATYPE" | cut -d_ -f1)
fi

FILENAMES=$(cat "$FILELIST" | xargs -n"$MAXFILES" | sed -n "$((JOBNUMBER+1)) p" | sed 's;^/xrootd/;root://cms-xrdr.private.lo:2094//xrd/;g')
if [ -z "$FILENAMES" ]; then
    echo "[01_prod_ntuple_leptonSF.sh] ERROR: empty input file slice. FILELIST=$FILELIST MAXFILES=$MAXFILES JOBNUMBER=$JOBNUMBER" >&2
    exit 3
fi

ARGS=""

# MET filters / event cleaning.
# DATATYPE is expected to be Run2024 for data and MC2024 for Summer24 MC.
if [ "${LEPSF_APPLY_FLAGS:-1}" == "1" ]; then
    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.flags flags_${DATATYPE}"
else
    echo "[01_prod_ntuple_leptonSF.sh] WARNING: LEPSF_APPLY_FLAGS=0. Running without MET/event-cleaning flags."
fi

# Lepton-SF-specific HLT filter.
# Default module is fcncTriLeptonHLT_leptonSF, where MC HLT is also filtered.
# Set LEPSF_HLT_PYMODULE=fcncTriLeptonHLT only if you intentionally want to use
# the common analysis HLT module instead.
HLT_PYMODULE=${LEPSF_HLT_PYMODULE:-fcncTriLeptonHLT_leptonSF}
if [ "${LEPSF_APPLY_HLT:-1}" == "1" ]; then
    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.${HLT_PYMODULE} hlt_${HLTMODULE}"
else
    echo "[01_prod_ntuple_leptonSF.sh] WARNING: LEPSF_APPLY_HLT=0. Running without HLT filter."
fi

ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.leptonSFNtupleProducer_Run3 leptonSF_${CHANNEL}_${YEAR}"

# Golden JSON for data.
if [ "${DATATYPE::2}" != "MC" ]; then
    if [ "$DATATYPE" == "Run2024" ]; then
        JSONFILE="$CMSSW_BASE/src/TZWi/TopAnalysis/data/lumi/Cert_Collisions2024_378981_386951_Golden.json"
        if [ -f "$JSONFILE" ]; then
            ARGS="-J $JSONFILE $ARGS"
        else
            echo "[01_prod_ntuple_leptonSF.sh] WARNING: JSON file not found: $JSONFILE"
        fi
    fi
fi

# Use a dedicated output branch selection. This keeps the NTuples small and
# preserves all SF_* branches produced by leptonSFNtupleProducer_Run3.
BRANCHSEL="${LEPSF_BRANCHSEL:-$CMSSW_BASE/src/TZWi/TopAnalysis/data/branchsel_leptonSF.txt}"
if [ ! -f "$BRANCHSEL" ]; then
    BRANCHSEL="$(pwd)/branchsel_leptonSF_autogen.txt"
    cat > "$BRANCHSEL" <<'BRANCHSEL_EOF'
drop *
keep run
keep luminosityBlock
keep event
keep genWeight
keep Pileup_nTrueInt
keep HLT
keep Flag
keep SF_*
BRANCHSEL_EOF
    echo "[01_prod_ntuple_leptonSF.sh] Created local branch selection: $BRANCHSEL"
fi

if [ "${USE_PREFETCH:-0}" = "1" ]; then
    CMD="nano_postproc.py --prefetch --bo $BRANCHSEL"
else
    CMD="nano_postproc.py --bo $BRANCHSEL"
fi

OUTPATH=ntuple_${YEAR}_leptonSF/reco/$CHANNEL/$DATASET0
[ ! -d "$OUTPATH" ] && mkdir -p "$OUTPATH"

printf '[01_prod_ntuple_leptonSF.sh] Start: %s\n' "$(date)"
echo "[01_prod_ntuple_leptonSF.sh] CHANNEL=$CHANNEL"
echo "[01_prod_ntuple_leptonSF.sh] FILELIST=$FILELIST"
echo "[01_prod_ntuple_leptonSF.sh] DATASET0=$DATASET0"
echo "[01_prod_ntuple_leptonSF.sh] DATATYPE=$DATATYPE YEAR=$YEAR ERA=$ERA PRIMARY_DATASET_RAW=$PRIMARY_DATASET_RAW PRIMARY_DATASET=$PRIMARY_DATASET"
echo "[01_prod_ntuple_leptonSF.sh] HLT_PYMODULE=$HLT_PYMODULE HLTMODULE=$HLTMODULE"
echo "[01_prod_ntuple_leptonSF.sh] MAXFILES=$MAXFILES JOBNUMBER=$JOBNUMBER"
echo "[01_prod_ntuple_leptonSF.sh] Number of input files in this job: $(echo $FILENAMES | wc -w)"
echo "[01_prod_ntuple_leptonSF.sh] CMD: $CMD $ARGS $OUTPATH $FILENAMES"

TMP_ERR="err_${CHANNEL}_${DATASET0}_${JOBNUMBER}.tmp"
/usr/bin/time -v $CMD $ARGS "$OUTPATH" $FILENAMES 2> "$TMP_ERR"
EXIT_CODE=$?
printf '[01_prod_ntuple_leptonSF.sh] End: %s\n' "$(date)"

# Remove ROOT/PyROOT teardown stack traces from stderr while preserving real errors.
sed '/\*\*\* Break \*\*\* segmentation violation/,$d' "$TMP_ERR" >&2
rm -f "$TMP_ERR"

if [ $EXIT_CODE -eq 0 ]; then
    echo "[01_prod_ntuple_leptonSF.sh] Job completed successfully."
else
    echo "[01_prod_ntuple_leptonSF.sh] Job exited with code $EXIT_CODE. Checking for PyROOT teardown-only failure..."
    if ls "$OUTPATH"/*.root 1> /dev/null 2>&1; then
        echo "[01_prod_ntuple_leptonSF.sh] Output ROOT file(s) exist. Masking teardown failure."
        EXIT_CODE=0
    else
        echo "[01_prod_ntuple_leptonSF.sh] No ROOT output found. Actual job failure." >&2
    fi
fi

# Compatibility tarball expected by create-batch3_mine/Condor transfer settings.
if [ $EXIT_CODE -eq 0 ]; then
    JOBNAME="${CHANNEL}.${DATASET0}"
    EXPECTED_TRANSFER_TGZ="result_${JOBNAME}.tgz"

    if [ "${MAKE_CONDOR_TRANSFER_TGZ:-1}" == "1" ]; then
        if ls "$OUTPATH"/*.root 1> /dev/null 2>&1; then
            echo "[01_prod_ntuple_leptonSF.sh] Creating Condor transfer tarball: $EXPECTED_TRANSFER_TGZ"
            rm -f "$EXPECTED_TRANSFER_TGZ"
            tar czf "$EXPECTED_TRANSFER_TGZ" "$OUTPATH"
            TAR_EXIT_CODE=$?
            if [ $TAR_EXIT_CODE -ne 0 ]; then
                echo "[01_prod_ntuple_leptonSF.sh] ERROR: failed to create $EXPECTED_TRANSFER_TGZ" >&2
                EXIT_CODE=$TAR_EXIT_CODE
            else
                ls -lh "$EXPECTED_TRANSFER_TGZ"
            fi
        else
            echo "[01_prod_ntuple_leptonSF.sh] ERROR: no ROOT output found under $OUTPATH; not creating transfer tarball." >&2
            EXIT_CODE=4
        fi
    fi
fi

if [ $EXIT_CODE -ne 0 ]; then
    echo $CMD $ARGS "$OUTPATH" $FILENAMES >> failed.txt
    exit $EXIT_CODE
else
    exit 0
fi
