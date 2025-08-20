#!/bin/bash

export BLUECRAB=/install/blue-crab-0.3.0/scripts/blue-crab
export SLOW5_DORADO=/install/slow5-dorado-0.8.3/bin/slow5-dorado
export DORADO_MODEL_PATH=/install/dorado-0.8.3/models/
export DORADO_MODEL=dna_r10.4.1_e8.2_400bps_sup@v5.0.0
export FSHD_PIPELINE_WRAPPER=/data/hasindu/realfshd/fshd_pipeline/fshd_pipeline_wrapper.sh

rm -f /data/hasindu/realfshd/test/data/positive/AS2603/realtime/*.log
rm -rf /data/hasindu/realfshd/test/data/positive/AS2603/realtime/sam /data/hasindu/realfshd/test/data/positive/AS2603/realtime/logs /data/hasindu/realfshd/test/data/positive/AS2603/realtime/slow5

scripts/realfreq.sh -n -m /data/hasindu/realfshd/test/data/positive/AS2603/realtime/
