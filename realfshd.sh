#!/bin/bash

REALPATH=$(dirname "$(readlink -f "$0")")

export BLUECRAB=/install/blue-crab-0.3.0/scripts/blue-crab
export SLOW5_DORADO=/install/slow5-dorado-0.8.3/bin/slow5-dorado
export DORADO_MODEL_PATH=/install/dorado-0.8.3/models/
export DORADO_MODEL=dna_r10.4.1_e8.2_400bps_sup@v5.0.0
export FSHD_PIPELINE_WRAPPER=${REALPATH}/fshd/fshd_pipeline_wrapper.sh

#DATA_DIR=${REALPATH}/test/data/positive/AS2603/realtime/
DATA_DIR=${REALPATH}/test/data/QGXXOX250638/
rm -f ${DATA_DIR}/*.log
rm -rf ${DATA_DIR}/sam ${DATA_DIR}/logs ${DATA_DIR}/slow5

${REALPATH}/scripts/realfreq.sh -n -m ${DATA_DIR}


