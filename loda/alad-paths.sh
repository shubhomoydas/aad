#!/bin/bash

if [[ "$OPR" == "1" ]]; then
    SCRIPT="alad.R"
elif [[ "$OPR" == "2" ]]; then
    SCRIPT="prep-alad-summary.R"
elif [[ "$OPR" == "3" ]]; then
    SCRIPT="ai2.R"
else
    SCRIPT="" # error
fi

SCRIPT_PATH=/scratch/code-submission/loda/R/$SCRIPT
FOLDER_PATH=/scratch/code-submission/datasets/anomaly/$DATASET/fullsamples
PLOTS_PATH=/scratch/code-submission/datasets/anomaly/$DATASET/fullplots
MODEL_PATH=/scratch/code-submission/datasets/anomaly/$DATASET/fullmodel
RESULTS_PATH=/scratch/code-submission/datasets/anomaly/$DATASET

if [[ "$OPR" == "3" ]]; then
    RESULTS_PATH="${RESULTS_PATH}/ai2results"
else
    RESULTS_PATH="${RESULTS_PATH}/fullresults"
fi

mkdir -p "${MODEL_PATH}"
mkdir -p "${RESULTS_PATH}"
mkdir -p "${PLOTS_PATH}"

