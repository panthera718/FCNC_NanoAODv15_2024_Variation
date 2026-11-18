#!/bin/bash

if [ $# != 4 ]; then
    echo "$#"
    echo "Usage: $0 FREl NanoAOD/2024/RD.Run2024/Run2024.JetMET/JetMET0.Run2024C-...NANOAOD.txt 10 0"
    echo "Usage: $0 FRMu NanoAOD/2024/RD.Run2024/Run2024.JetMET/JetMET0.Run2024C-...NANOAOD.txt 10 0"
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

if [ "$CHANNEL" != "FRMu" ] && [ "$CHANNEL" != "FREl" ]; then
    echo "[01_prod_ntuple_FR.sh] ERROR: supported channels are FRMu and FREl, got CHANNEL=$CHANNEL"
    exit 2
fi

DATASET0=`basename $FILELIST | sed -e 's;.txt;;g'`
DATASET='/'`echo $DATASET0 | sed -e 's;\.;/;g'`
ERA=$(echo $DATASET0 | cut -d. -f2 | cut -d- -f1 | sed -e 's;NanoAOD;;g')
DATATYPE=$(basename $(dirname $FILELIST) | cut -d. -f1)
YEAR=${DATATYPE:(-4)}
PRIMARY_DATASET=$(echo $DATASET | cut -d/ -f2 | sed -E 's/[01]$//')

# Guard against accidental stream mixing. Nominal FR must run on JetMET data.
# EGamma/Muon may be allowed only for explicit cross-checks with ALLOW_LEPTON_FR_STREAMS=1.
if [ ${DATATYPE::3} == "Run" ]; then
    if [ "$PRIMARY_DATASET" != "JetMET" ]; then
        if [ "${ALLOW_LEPTON_FR_STREAMS:-0}" != "1" ]; then
            echo "[01_prod_ntuple_FR.sh] ERROR: nominal FR must run on JetMET0/JetMET1, got $PRIMARY_DATASET from $DATASET0"
            echo "[01_prod_ntuple_FR.sh]        Set ALLOW_LEPTON_FR_STREAMS=1 only for non-nominal EGamma/Muon cross-checks."
            exit 2
        fi
    fi
fi

FILENAMES=$(cat $FILELIST | xargs -n$MAXFILES | sed -n "$((JOBNUMBER+1)) p" | sed 's;^/xrootd/;root://cms-xrdr.private.lo:2094//xrd/;g')
if [ -z "$FILENAMES" ]; then
    echo "[01_prod_ntuple_FR.sh] ERROR: no input files selected for JOBNUMBER=$JOBNUMBER FILELIST=$FILELIST"
    exit 3
fi

BRANCHSEL="$CMSSW_BASE/src/TZWi/TopAnalysis/data/branchsel_FR.txt"
if [ ! -f "$BRANCHSEL" ]; then
    echo "[01_prod_ntuple_FR.sh] WARNING: $BRANCHSEL does not exist. Running without --bo branch selection."
    CMD="nano_postproc.py --prefetch"
else
    CMD="nano_postproc.py --prefetch --bo $BRANCHSEL"
fi

# Nominal fake-rate trigger configuration. The Python producer reads these env vars.
export FR_TRIGGER_YAML=${FR_TRIGGER_YAML:-$CMSSW_BASE/src/TZWi/TopAnalysis/data/HLT/2024_FR_JetMET_Trigger.yaml}
export FR_TRIGGER_KEY=${FR_TRIGGER_KEY:-Run2024.JetMET_FR}
export FR_MIN_LEAD_JET_PT=${FR_MIN_LEAD_JET_PT:-250.0}
export FR_MAX_MET_PLUS_MT=${FR_MAX_MET_PLUS_MT:-50.0}
export FR_MAX_LOOSE_LEPTONS=${FR_MAX_LOOSE_LEPTONS:-1}

# Analysis-consistent jet and UParTAK4 b-tag counting stored in the FR NTuple.
# The trigger jet remains pT>FR_MIN_LEAD_JET_PT; it is not used as the b-jet pT cut.
export FR_ANALYSIS_MIN_JET_PT=${FR_ANALYSIS_MIN_JET_PT:-30.0}
export FR_ANALYSIS_MAX_JET_ABS_ETA=${FR_ANALYSIS_MAX_JET_ABS_ETA:-4.7}
export FR_BTAG_MIN_JET_PT=${FR_BTAG_MIN_JET_PT:-20.0}
export FR_BTAG_MAX_JET_ABS_ETA=${FR_BTAG_MAX_JET_ABS_ETA:-2.5}
export FR_JET_LEPTON_DR=${FR_JET_LEPTON_DR:-0.4}
export FR_BTAG_WP=${FR_BTAG_WP:-L}
# Leave FR_BTAG_WP_VALUE empty to use the official 2024 value for FR_BTAG_WP.
export FR_BTAG_WP_VALUE=${FR_BTAG_WP_VALUE:-}
export FR_APPLY_OFFICIAL_JET_ID=${FR_APPLY_OFFICIAL_JET_ID:-1}
export FR_JETID_WP=${FR_JETID_WP:-AK4PUPPI_Tight}
export FR_JETID_JSON=${FR_JETID_JSON:-$CMSSW_BASE/src/TZWi/TopAnalysis/data/JetID/2024/jetid.json.gz}

OUTPATH="ntuple_${YEAR}_FR/reco/${CHANNEL}/${DATASET0}"
mkdir -p "$OUTPATH"

# The producer applies event-cleaning filters, JetMET FR trigger, one-probe-lepton
# selection, and the promptMVA numerator definition internally.
ARGS="-I TZWi.TopAnalysis.postprocessing.promptMVAFakeRateNtupleProducer_Run3 FRProducer_${CHANNEL}_${YEAR}"

# Debug-only: disable internal HLT. Do not use this for nominal fake-rate maps.
if [ "${FR_REQUIRE_TRIGGER:-1}" == "0" ]; then
    ARGS="-I TZWi.TopAnalysis.postprocessing.promptMVAFakeRateNtupleProducer_Run3 FRProducerNoHLT_${CHANNEL}_${YEAR}"
fi

echo "[01_prod_ntuple_FR.sh] CHANNEL=$CHANNEL"
echo "[01_prod_ntuple_FR.sh] DATASET0=$DATASET0"
echo "[01_prod_ntuple_FR.sh] DATATYPE=$DATATYPE YEAR=$YEAR ERA=$ERA PRIMARY_DATASET=$PRIMARY_DATASET"
echo "[01_prod_ntuple_FR.sh] OUTPATH=$OUTPATH"
echo "[01_prod_ntuple_FR.sh] FR_TRIGGER_YAML=$FR_TRIGGER_YAML"
echo "[01_prod_ntuple_FR.sh] FR_TRIGGER_KEY=$FR_TRIGGER_KEY"
echo "[01_prod_ntuple_FR.sh] FR_MIN_LEAD_JET_PT=$FR_MIN_LEAD_JET_PT"
echo "[01_prod_ntuple_FR.sh] FR_MAX_MET_PLUS_MT=$FR_MAX_MET_PLUS_MT"
echo "[01_prod_ntuple_FR.sh] FR_ANALYSIS_MIN_JET_PT=$FR_ANALYSIS_MIN_JET_PT"
echo "[01_prod_ntuple_FR.sh] FR_ANALYSIS_MAX_JET_ABS_ETA=$FR_ANALYSIS_MAX_JET_ABS_ETA"
echo "[01_prod_ntuple_FR.sh] FR_BTAG_MIN_JET_PT=$FR_BTAG_MIN_JET_PT"
echo "[01_prod_ntuple_FR.sh] FR_BTAG_MAX_JET_ABS_ETA=$FR_BTAG_MAX_JET_ABS_ETA"
echo "[01_prod_ntuple_FR.sh] FR_JET_LEPTON_DR=$FR_JET_LEPTON_DR"
echo "[01_prod_ntuple_FR.sh] FR_BTAG_WP=$FR_BTAG_WP"
echo "[01_prod_ntuple_FR.sh] FR_BTAG_WP_VALUE=${FR_BTAG_WP_VALUE:-<official>}"
echo "[01_prod_ntuple_FR.sh] FR_JETID_WP=$FR_JETID_WP"
echo "[01_prod_ntuple_FR.sh] FR_JETID_JSON=$FR_JETID_JSON"
echo "[01_prod_ntuple_FR.sh] FILES=$FILENAMES"
echo "[01_prod_ntuple_FR.sh] ARGS=$ARGS"

$CMD "$OUTPATH" $FILENAMES $ARGS
RET=$?
if [ $RET -ne 0 ]; then
    echo "[01_prod_ntuple_FR.sh] ERROR: nano_postproc.py failed with exit code $RET"
    exit $RET
fi

# create-batch3_mine expects a result tarball in the running directory.
tar czf "result_${CHANNEL}_${DATASET0}_${JOBNUMBER}.tgz" "$OUTPATH"
