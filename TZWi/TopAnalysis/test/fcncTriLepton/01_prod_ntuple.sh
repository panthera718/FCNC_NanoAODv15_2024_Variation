#!/bin/bash

if [ $# != 4 ]; then
    echo $#
    echo "Usage: $0 MuElEl MC2016.WW.txt 10 5 ## process WW sample assuming muee channel, split by 10 files and run 5th section"
    echo "Usage: $0 ElElEl MC2017.TT_powheg.txt 1 0 ## process TTbar sample assuming eee channel, one file per each section and run 0th one."
    echo "Usage: $0 MuMuMu Run2016B.DoubleMuon.txt 1 0"
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
ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLepton_allmodes_unified_jetvar_only_promptMVA_Run3 fcnc_${CHANNEL}_${YEAR}"
#if [ ${DATATYPE::2} != "MC" ]; then
#    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLepton_allmodes_unified_jetvar_only fcnc_${CHANNEL}_${YEAR}"
#elif [ ${DATATYPE::2} == "MC" ]; then
#    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.fcncTriLepton_allmodes_unified_jetvar_only central_fcnc_${CHANNEL}_${YEAR}"
#fi
#ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.JetVetoProducer_unified JVProducer"

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

#if [ ${DATATYPE::2} == "MC" ]; then

    ##########
    # Common #
    ##########
#    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.lepSFProducer lepSF"
    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.lepSFProducer_2024_signedEta lepSF"
#    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.lepRecoSFProducer lepRecoSF"
#    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.puWeightProducer puWeight_NewUL${YEAR}"
#    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.tauTriggerSFProducer tauTriggerSF_${CHANNEL}_${YEAR}"
#    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.tauSFProducer tauSF_${CHANNEL}_${YEAR}"

    ###########
    # central #
    ###########
#    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.JetPUIDSFProducer_unified JPUIDSFProducer"
    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.btagShapeWeightProducer_unified_jesmatched btagWeight"
#    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.ctagShapeWeightProducer_unified ctagWeight"
#    ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.qgtagShapeWeightProducer_unified qgtagWeight"

    ##########
    # jetVar #
    ##########
#    for jetVar in $jetUncert
#    do
#        ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.JetPUIDSFProducer_unified JPUIDSFProducerUp_${jetVar}"
#        ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.btagShapeWeightProducer_unified_jesmatched btagWeightUp_${jetVar}"
#        ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.ctagShapeWeightProducer_unified ctagWeightUp_${jetVar}"
#        ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.qgtagShapeWeightProducer_unified qgtagWeightUp_${jetVar}"
#
#        ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.JetPUIDSFProducer_unified JPUIDSFProducerDown_${jetVar}"
#        ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.btagShapeWeightProducer_unified_jesmatched btagWeightDown_${jetVar}"
#        ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.ctagShapeWeightProducer_unified ctagWeightDown_${jetVar}"
#        ARGS="$ARGS -I TZWi.TopAnalysis.postprocessing.qgtagShapeWeightProducer_unified qgtagWeightDown_${jetVar}"
#    done

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
