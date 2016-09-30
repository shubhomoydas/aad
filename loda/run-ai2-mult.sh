#!/bin/bash

# To run:
# bash ./run-ai2-mult.sh <dataset> <reruns> <tau>
#
# Examples:
# bash ./run-ai2-mult.sh toy 10 0.03

DATASET=$1
RERUNS=$2
TAU=$3

OPR=3 # OPR: 1(alad) | 2(alad-summary) | 3(ai2)
source ./alad-paths.sh
echo "Rscript $SCRIPT_PATH"

REPS=1
MAX_BUDGET=100
BUDGET=0
TOPK=0

if [[ "$DATASET" == "yeast" ]]; then
    BUDGET=60 # just to make this similar to other larger datasets
fi

echo "Options: ai2 reps-${REPS} budget-$BUDGET tau-${TAU} topK-$TOPK"

Rscript $SCRIPT_PATH --dataset=$DATASET --reps=$REPS --reruns=$RERUNS --budget=$BUDGET --maxbudget=$MAX_BUDGET --topK=$TOPK --tau=$TAU --filedir=$FOLDER_PATH --cachedir=$MODEL_PATH --resultsdir=$RESULTS_PATH --plotsdir=$PLOTS_PATH
