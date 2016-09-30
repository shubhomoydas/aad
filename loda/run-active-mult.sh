#!/bin/bash

# To run:
# bash ./run-active.sh <dataset> <inferencetype[1|2|3]> <unif_prior[0|1]> <tau> <querytype[1|3]>
#
# inferencetype: 1-Simple, 2-Simple_Optim(same as simple), 3-AATP with pairwise constraints
# querytype: 1-Top, 3-quantile
#
# IMPORTANT: Currently *only* inferencetype 3 and querytype 1 are used
#
# Examples:
# bash ./run-active-mult.sh toy 3 1 0.03 1

DATASET=$1
INFERENCE_TYPE=$2
UNIF_PRIOR_IND=$3
TAU=$4
Q_TYPE=$5

SIGMA2=0.5

OPR=1 # OPR: 1(alad) | 2(alad-summary)
source ./alad-paths.sh
echo "Rscript $SCRIPT_PATH"

PRIOR_INDS=( "1" )
MEAN_REL_LOSS_INDS=( "0" )
QUERY_TYPES=( "${Q_TYPE}" )
BATCH_INDS=( "0" )
PSEUDO_ALWAYS=0

REPS=1
RERUNS=10
BUDGET=0
MAX_BUDGET=100
TOPK=0

if [[ "$PSEUDO_ALWAYS" == "1" ]]; then
    PSEUDO_OPT="--pseudoanomrank_always"
else
    PSEUDO_OPT=""
fi

for WITH_PRIOR_IND in "${PRIOR_INDS[@]}"
do
    for MEAN_REL_LOSS_IND in "${MEAN_REL_LOSS_INDS[@]}"
    do
        for QUERY_TYPE in "${QUERY_TYPES[@]}"
        do
            for BATCH_IND in "${BATCH_INDS[@]}"
            do

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

            echo "Options: inference-${INFERENCE_TYPE} ${WITH_PRIOR} query-${QUERY_TYPE} ${MEAN_REL_LOSS} ${BATCH} reps-${REPS} budget-$BUDGET maxbudget-${MAX_BUDGET} tau-${TAU} topK-$TOPK ${UNIF_PRIOR} ${PSEUDO_OPT}"

            Rscript $SCRIPT_PATH --dataset=$DATASET --inferencetype=$INFERENCE_TYPE --sigma2=$SIGMA2 $WITH_PRIOR $UNIF_PRIOR $MEAN_REL_LOSS $BATCH --reps=$REPS --reruns=$RERUNS --budget=$BUDGET --maxbudget=${MAX_BUDGET} --topK=$TOPK --tau=$TAU --filedir=$FOLDER_PATH --querytype=$QUERY_TYPE ${PSEUDO_OPT} --cachedir=$MODEL_PATH --resultsdir=$RESULTS_PATH --plotsdir=$PLOTS_PATH

            done # BATCH_IND
        done # QUERY_TYPE
    done # MEAN_REL_IND
done # PRIOR_IND
