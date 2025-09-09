#!/bin/bash

REALPATH=$(dirname "$(readlink -f "$0")")
${REALPATH}/../fshd/fshd_pipeline_wrapper.sh ${REALPATH}/data/PGXXXX230339_reads_20k.bam SAMPLE
