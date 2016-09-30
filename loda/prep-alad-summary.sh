#!/bin/bash

# To run:
# bash ./prep-alad-summary.sh <dataset> <inferencetype[1|2|3]> <reps[1-20]> <tau>
#
# inferencetype: 1-Simple, 2-Simple_Optim(same as simple), 3-AATP with pairwise constraints
#
# Examples:
# bash ./prep-alad-summary.sh toy 3 1 0.03

DATASET=$1
INFERENCE_TYPE=$2
REPS=$3
TAU=$4

SIGMA2=0.5

OPR=2 # OPR: 1(alad) | 2(alad-summary)
source ./alad-paths.sh
echo "Rscript $SCRIPT_PATH"

RERUNS=10
BUDGET=0
MAX_BUDGET=100
TOPK=0
TAU=0.03

WITH_PRIOR_IND=1
MEAN_REL_LOSS_IND=0
QUERY_TYPE=1
BATCH_IND=0
UNIF_PRIOR_IND=1
PSEUDO_ALWAYS=0

if [[ "$PSEUDO_ALWAYS" == "1" ]]; then
    PSEUDO_OPT="--pseudoanomrank_always"
else
    PSEUDO_OPT=""
fi

if [[ "$WITH_PRIOR_IND" == "1" ]]; then
    WITH_PRIOR="--withprior"
else
    WITH_PRIOR=""
fi

if [[ "$MEAN_REL_LOSS_IND" == "1" ]]; then
    MEAN_REL_LOSS="--withmeanrelativeloss"
else
    MEAN_REL_LOSS=""
fi

if [[ "$BATCH_IND" == "1" ]]; then
    BATCH="--batch"
else
    BATCH=""
fi

if [[ "$UNIF_PRIOR_IND" == "1" ]]; then
    UNIF_PRIOR="--unifprior"
else
    UNIF_PRIOR=""
fi

echo "Options: inference-${INFERENCE_TYPE} query-${QUERY_TYPE} reps-${REPS} tau-$TAU budget-$BUDGET"

Rscript $SCRIPT_PATH --dataset=$DATASET --inferencetype=$INFERENCE_TYPE --sigma2=$SIGMA2 $WITH_PRIOR $UNIF_PRIOR $MEAN_REL_LOSS $BATCH --reps=$REPS --reruns=$RERUNS --budget=$BUDGET --maxbudget=$MAX_BUDGET --topK=$TOPK --tau=$TAU --filedir=$FOLDER_PATH --querytype=$QUERY_TYPE ${PSEUDO_OPT} --cachedir=$MODEL_PATH --resultsdir=$RESULTS_PATH --plotsdir=$PLOTS_PATH
