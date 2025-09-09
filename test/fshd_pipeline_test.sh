#!/bin/bash

REALPATH=$(dirname "$(readlink -f "$0")")
${REALPATH}/../fshd/fshd_pipeline_wrapper.sh ${REALPATH}/data/positive/AS2603/QGXHXX240243_supmeth_unaligned.bam SAMPLE
