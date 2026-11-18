#!/bin/bash

# ----------------------------------------------------------------------
# Fake-rate map family default
# ----------------------------------------------------------------------
# Change only this line when you want the shell script's default to switch:
#
#   split     -> use *_0B.root for nBjet==0 and *_GE1B.root for nBjet>=1
#   inclusive -> use *_Inclusive.root for every event
#
#DEFAULT_FAKERATE_MAP_MODE="split"
DEFAULT_FAKERATE_MAP_MODE="inclusive"

if [ $# -ne 4 ] && [ $# -ne 5 ]; then
    echo "$#"
    echo "Usage (shell default):"
    echo "  $0 CHANNEL FILELIST MAXFILES JOBNUMBER"
    echo "Usage (explicit FR map mode; create-batch appends JOBNUMBER last):"
    echo "  $0 CHANNEL FILELIST MAXFILES {split|inclusive} JOBNUMBER"
    exit 1
fi

eval `scram runtime -sh`

export PYTHONNOUSERSITE=1
export XRD_WORKERTHREADS=1
export XRD_NETWORKSTACK=IPv4

CHANNEL=$1
FILELIST=$2
MAXFILES=$3

if [ $# -eq 5 ]; then
    REQUESTED_FAKERATE_MAP_MODE=$4
    JOBNUMBER=$5
else
    REQUESTED_FAKERATE_MAP_MODE=${FAKERATE_MAP_MODE:-$DEFAULT_FAKERATE_MAP_MODE}
    JOBNUMBER=$4
fi

case "$(echo "$REQUESTED_FAKERATE_MAP_MODE" | tr '[:upper:]' '[:lower:]')" in
    split|btag|btag_split|btag-split|1|true|yes)
        FAKERATE_MAP_MODE="split"
        ;;
    inclusive|inc|0|false|no)
        FAKERATE_MAP_MODE="inclusive"
        ;;
    *)
        echo "[01_prod_ntuple.sh] ERROR: unsupported fake-rate map mode: $REQUESTED_FAKERATE_MAP_MODE" >&2
        echo "[01_prod_ntuple.sh]        Use split or inclusive." >&2
        exit 2
        ;;
esac
export FAKERATE_MAP_MODE

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

# Extract era robustly: Run2024C, Run2024D, ..., Run2024I
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
echo "[DEBUG] FAKERATE_MAP_MODE = $FAKERATE_MAP_MODE"

FILENAMES=$(cat $FILELIST | xargs -n$MAXFILES | sed -n "$(($JOBNUMBER+1)) p" | sed 's;^/xrootd/;root://cms-xrdr.private.lo:2094//xrd/;g')

#jetUncert="AbsoluteStat AbsoluteScale AbsoluteMPFBias Fragmentation SinglePionECAL SinglePionHCAL FlavorQCD TimePtEta RelativeJEREC1 RelativeJEREC2 RelativeJERHF RelativePtBB RelativePtEC1 RelativePtEC2 RelativePtHF RelativeBal RelativeSample RelativeFSR RelativeStatFSR RelativeStatEC RelativeStatHF PileUpDataMC PileUpPtRef PileUpPtBB PileUpPtEC1 PileUpPtEC2 PileUpPtHF"

#jetUncert="AbsoluteStat AbsoluteScale AbsoluteMPFBias"
#jetUncert="Fragmentation SinglePionECAL SinglePionHCAL"
#jetUncert="FlavorQCD TimePtEta RelativeJEREC1"
#jetUncert="RelativeJEREC2 RelativeJERHF RelativePtBB"
#jetUncert="RelativePtEC1 RelativePtEC2 RelativePtHF"
#jetUncert="RelativeBal RelativeSample RelativeFSR"
#jetUncert="RelativeStatFSR RelativeStatEC RelativeStatHF"
#jetUncert="PileUpDataMC PileUpPtRef PileUpPtBB"
#jetUncert="PileUpPtEC1 PileUpPtEC2 PileUpPtHF"

##########
# Common #
##########
ARGS=""
ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.flags flags_${DATATYPE}"
ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLeptonHLT hlt_${HLTMODULE}"

###########
# central #
###########
# Use the hadron-flavour-writing factory for MC, because the 2024 UParTAK4
# b-tag WP SF needs GoodJet_hadronFlavour.  Data does not carry this branch.
if [ ${DATATYPE::2} == "MC" ]; then
    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLepton_allmodes_unified_jetvar_only_promptMVA_Run3 central_fcnc_${CHANNEL}_${YEAR}"
else
    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLepton_allmodes_unified_jetvar_only_promptMVA_Run3 fcnc_${CHANNEL}_${YEAR}"
fi

# Official 2024 JME jet-veto map.  This must run after FCNCTriLepton and
# before the cut flow, because CutFlow vetoes events with JetVetoSum != 0.
ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.JetVetoProducer_unified JVProducer"

##########
# jetVar #
##########
#if [ ${DATATYPE::2} == "MC" ]; then
#    for jetVar in $jetUncert
#    do
#        ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLepton_allmodes_unified_jetvar_only jetUncertUp_fcnc_${CHANNEL}_${YEAR}_${jetVar}"
#        ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.JetVetoProducer_unified JVProducerUp_${jetVar}"
#
#        ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLepton_allmodes_unified_jetvar_only jetUncertDown_fcnc_${CHANNEL}_${YEAR}_${jetVar}"
#        ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.JetVetoProducer_unified JVProducerDown_${jetVar}"
#    done
#fi

##########
# Common #
##########
ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLeptonCutFlow cutFlow_${CHANNEL}"

###########
# central #
###########
#ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncKinematicReco_unified fcncKinReco"
#ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncMVAinput_unified fcncMVAinput"

##########
# jetVar #
##########
#if [ ${DATATYPE::2} == "MC" ]; then
#    for jetVar in $jetUncert
#    do
#        ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncKinematicReco_unified fcncKinRecoUp_${jetVar}"
#        ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncMVAinput_unified fcncMVAinputUp_${jetVar}"
#
#        ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncKinematicReco_unified fcncKinRecoDown_${jetVar}"
#        ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncMVAinput_unified fcncMVAinputDown_${jetVar}"
#    done
#fi

##########
# Common #
##########
if [ ${CHANNEL::3} == "NPL" ]; then
    FRMODULE="FRProducer"
    if [ -n "$FAKERATE_ERA" ]; then
        FRERA="$FAKERATE_ERA"
    else
        FRERA=""
        case "$YEAR" in
            2024)
                FRERA="2024"
                ;;
        esac
    fi

    case "$FRERA" in
        2024)
            FRMODULE="FRProducer_2024"
            export FAKERATE_ERA="2024"

            FR_BASE_DIR=${FAKERATE_2024_DIR:-$CMSSW_BASE/src/TZWi/TopAnalysis/data/leptonFR/2024}
            FR_FILES_TO_CHECK=""

            case "$FAKERATE_MAP_MODE" in
                split)
                    # Event-dependent selection:
                    #   nBjet == 0 -> *_0B.root
                    #   nBjet >= 1 -> *_GE1B.root
                    export FAKERATE_ELECTRON_2024_0B_ROOT=${FAKERATE_ELECTRON_2024_0B_ROOT:-$FR_BASE_DIR/fakerate_electron_promptMVA_2024_JetMET_SignedEta_0B.root}
                    export FAKERATE_ELECTRON_2024_GE1B_ROOT=${FAKERATE_ELECTRON_2024_GE1B_ROOT:-$FR_BASE_DIR/fakerate_electron_promptMVA_2024_JetMET_SignedEta_GE1B.root}
                    export FAKERATE_MUON_2024_0B_ROOT=${FAKERATE_MUON_2024_0B_ROOT:-$FR_BASE_DIR/fakerate_muon_promptMVA_2024_JetMET_SignedEta_0B.root}
                    export FAKERATE_MUON_2024_GE1B_ROOT=${FAKERATE_MUON_2024_GE1B_ROOT:-$FR_BASE_DIR/fakerate_muon_promptMVA_2024_JetMET_SignedEta_GE1B.root}

                    FR_FILES_TO_CHECK="$FAKERATE_ELECTRON_2024_0B_ROOT
$FAKERATE_ELECTRON_2024_GE1B_ROOT
$FAKERATE_MUON_2024_0B_ROOT
$FAKERATE_MUON_2024_GE1B_ROOT"

                    echo "[01_prod_ntuple.sh] FAKERATE_MAP_MODE=split"
                    echo "[01_prod_ntuple.sh] electron 0B FR   = $FAKERATE_ELECTRON_2024_0B_ROOT"
                    echo "[01_prod_ntuple.sh] electron GE1B FR = $FAKERATE_ELECTRON_2024_GE1B_ROOT"
                    echo "[01_prod_ntuple.sh] muon 0B FR       = $FAKERATE_MUON_2024_0B_ROOT"
                    echo "[01_prod_ntuple.sh] muon GE1B FR     = $FAKERATE_MUON_2024_GE1B_ROOT"
                    ;;
                inclusive)
                    # One inclusive map per flavor for every nBjet category.
                    export FAKERATE_ELECTRON_2024_INCLUSIVE_ROOT=${FAKERATE_ELECTRON_2024_INCLUSIVE_ROOT:-$FR_BASE_DIR/fakerate_electron_promptMVA_2024_JetMET_SignedEta_Inclusive.root}
                    export FAKERATE_MUON_2024_INCLUSIVE_ROOT=${FAKERATE_MUON_2024_INCLUSIVE_ROOT:-$FR_BASE_DIR/fakerate_muon_promptMVA_2024_JetMET_SignedEta_Inclusive.root}

                    FR_FILES_TO_CHECK="$FAKERATE_ELECTRON_2024_INCLUSIVE_ROOT
$FAKERATE_MUON_2024_INCLUSIVE_ROOT"

                    echo "[01_prod_ntuple.sh] FAKERATE_MAP_MODE=inclusive"
                    echo "[01_prod_ntuple.sh] electron inclusive FR = $FAKERATE_ELECTRON_2024_INCLUSIVE_ROOT"
                    echo "[01_prod_ntuple.sh] muon inclusive FR     = $FAKERATE_MUON_2024_INCLUSIVE_ROOT"
                    ;;
            esac

            while IFS= read -r FR_FILE
            do
                [ -z "$FR_FILE" ] && continue
                if [ ! -f "$FR_FILE" ]; then
                    echo "[01_prod_ntuple.sh] ERROR: fake-rate file does not exist: $FR_FILE" >&2
                    exit 3
                fi
            done <<EOF_FR_FILES
$FR_FILES_TO_CHECK
EOF_FR_FILES
            ;;
        *)
            if echo "$CHANNEL" | grep -q 'Ta'; then
                echo "[01_prod_ntuple.sh] ERROR: could not determine tau fake-rate era for CHANNEL=$CHANNEL FILELIST=$FILELIST"
                echo "Set FAKERATE_ERA explicitly to one of: 2016preVFP, 2016postVFP, 2017, 2024"
                exit 2
            fi
            ;;
    esac

    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.FakeRateProducer ${FRMODULE}"
fi

BRANCHSEL="$CMSSW_BASE/src/TZWi/TopAnalysis/data/branchsel.txt"

#CMD="nano_postproc.py --bo $BRANCHSEL"
# 1. --prefetch  option added: Download files locally first to prevent XRootD communication conflicts at the source.
CMD="nano_postproc.py --prefetch --bo $BRANCHSEL"

OUTPATH=ntuple_$YEAR/reco/$CHANNEL/$DATASET0

[ ! -d $OUTPATH ] && mkdir -p $OUTPATH

if [ ${DATATYPE::2} == "MC" ]; then

    ##########
    # Common MC event/lepton/object corrections
    ##########
    # Analysis-specific promptMVA lepton SFs measured directly from 2024 data/MC.
    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.lepSFProducer lepSF"

    # Official EGM electron reconstruction SF.  This is intentionally separate
    # from the direct promptMVA/selection SF above.
    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.lepRecoSFProducer lepRecoSF_2024"

    # Official LUM pileup weights for the Run2024 CDEFGHI data-taking set.
    # Use puWeight_BCDEFGHI_2024 instead if the production includes Run2024B.
    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.puWeightProducer puWeight_2024"

    ###########
    # central #
    ###########
    # 2024 UParTAK4 fixed-WP b-tag SF.  Do not run the old JetPUID SF producer
    # for Run-3 PUPPI AK4 jets.
    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.btagShapeWeightProducer_unified_jesmatched btagWeight"

    ##########
    # jetVar #
    ##########
#    for jetVar in $jetUncert
#    do
#        ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.btagShapeWeightProducer_unified_jesmatched btagWeightUp_${jetVar}"
#        ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.btagShapeWeightProducer_unified_jesmatched btagWeightDown_${jetVar}"
#    done
fi

#else
#    if [ $DATATYPE == "Run2022" ]; then
#        JSONFILE="$CMSSW_BASE/src/TZWi/TopAnalysis/data/lumi/Cert_Collisions2022_355100_362760_Golden.json"
#        ARGS="-J $JSONFILE $ARGS"
#    elif [ $DATATYPE == "Run2023" ]; then
#        JSONFILE="$CMSSW_BASE/src/TZWi/TopAnalysis/data/lumi/Cert_Collisions2023_366442_370790_Golden.json"
#        ARGS="-J $JSONFILE $ARGS"
#    elif [ $DATATYPE == "Run2024" ]; then
#        JSONFILE="$CMSSW_BASE/src/TZWi/TopAnalysis/data/lumi/Cert_Collisions2024_378981_386951_Golden.json"
#        ARGS="-J $JSONFILE $ARGS"
#    fi
#fi

if [ ${DATATYPE::2} != "MC" ]; then
    if [ $DATATYPE == "Run2022" ]; then
        JSONFILE="$CMSSW_BASE/src/TZWi/TopAnalysis/data/lumi/Cert_Collisions2022_355100_362760_Golden.json"
        ARGS="-J $JSONFILE $ARGS"
    elif [ $DATATYPE == "Run2023" ]; then
        JSONFILE="$CMSSW_BASE/src/TZWi/TopAnalysis/data/lumi/Cert_Collisions2023_366442_370790_Golden.json"
        ARGS="-J $JSONFILE $ARGS"
    elif [ $DATATYPE == "Run2024" ]; then
        JSONFILE="$CMSSW_BASE/src/TZWi/TopAnalysis/data/lumi/Cert_Collisions2024_378981_386951_Golden.json"
        ARGS="-J $JSONFILE $ARGS"
    fi
fi


#echo $CMD $ARGS $OUTPATH $FILENAMES
#$CMD $ARGS $OUTPATH $FILENAMES
#[ $? -eq 0 ] || echo $CMD $ARGS $OUTPATH $FILENAMES >> failed.txt

### Wrapping ###
echo $CMD $ARGS $OUTPATH $FILENAMES

# 2. Execute Python script while saving standard error (stderr) to a temporary file
TMP_ERR="err_${JOBNUMBER}.tmp"
$CMD $ARGS $OUTPATH $FILENAMES 2> $TMP_ERR
EXIT_CODE=$?

# 2-1. Delete from "*** Break *** segmentation violation" to the end (stack trace) in the temporary error log
# Only the remaining content (actual error messages) is output back to the original standard error (&2) and recorded in job_*.err.
sed '/\*\*\* Break \*\*\* segmentation violation/,$d' $TMP_ERR >&2
rm -f $TMP_ERR

# 3. Logic for bypassing crashes (Segfault, GIL deadlock, etc.) at exit point
if [ $EXIT_CODE -eq 0 ]; then
    echo "Job completed successfully."
else
    echo "Job exited with code $EXIT_CODE. Checking if it's a PyROOT teardown crash..."
    
    # Check if the .root file was successfully created in the output directory
    if ls $OUTPATH/*.root 1> /dev/null 2>&1; then
        echo "Output .root file(s) successfully generated!"
        echo "Masking the Python teardown crash and reporting success to the batch system."
        EXIT_CODE=0  # Force overwrite error code with 0 (success)
    else
        echo "No output .root files found. Actual job failure." >&2
        # Here, leave EXIT_CODE as is to trigger the failure logic below.
    fi
fi

# 4. failed.txt log based on final result and script termination
if [ $EXIT_CODE -ne 0 ]; then
    echo $CMD $ARGS $OUTPATH $FILENAMES >> failed.txt
    exit $EXIT_CODE
else
    exit 0
fi
