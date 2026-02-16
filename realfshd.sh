#!/bin/bash

REALPATH=$(dirname "$(readlink -f "$0")")

export BLUECRAB=${REALPATH}/fshd/install/blue-crab-v0.4.0/blue-crab
export SLOW5_DORADO=${REALPATH}/fshd/install/slow5-dorado-0.8.3/bin/slow5-dorado
export DORADO_MODEL_PATH=${REALPATH}/fshd/install/slow5-dorado-0.8.3/models/
export DORADO_MODEL_FAST=dna_r10.4.1_e8.2_400bps_fast@v5.0.0
export DORADO_MODEL_SUP=dna_r10.4.1_e8.2_400bps_sup@v5.0.0
export FSHD_PRE=${REALPATH}/fshd/pre_mapping_n_rebasecall.sh
export FSHD_PIPELINE_WRAPPER=${REALPATH}/fshd/fshd_pipeline_wrapper.sh
export FSHD_POST=${REALPATH}/fshd/fshd_post_merge.sh

export FSHDTOOLS=${REALPATH}/fshdtools

export CUDA_VISIBLE_DEVICES=1

export PATH=$PATH://install/minimap2-2.26/:/install/samtools-1.18/:/install/slow5tools-v1.3.0/


${REALPATH}/scripts/realfreq.sh $*

