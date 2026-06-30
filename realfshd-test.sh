#!/bin/bash

REALPATH=$(dirname "$(readlink -f "$0")")

: "${BLUECRAB:=$HOME/tools/venvs/bluecrab/bin/blue-crab}"
export BLUECRAB
export SLOW5_DORADO=/install/slow5-dorado-1.1.1/bin/slow5-dorado
export DORADO_MODEL_PATH=/install/slow5-dorado-1.1.1/models/
export DORADO_MODEL_FAST=dna_r10.4.1_e8.2_400bps_fast@v5.2.0
export DORADO_MODEL_SUP=dna_r10.4.1_e8.2_400bps_sup@v5.2.0
export FSHD_PRE="${REALPATH}/fshd/pre_mapping_n_rebasecall.sh"
export FSHD_PIPELINE_WRAPPER="${REALPATH}/fshd/fshd_pipeline_wrapper.sh"
export FSHD_POST="${REALPATH}/fshd/fshd_post_merge.sh"
export FSHDTOOLS="${REALPATH}/fshdtools"

# Set maximum performance mode and lock all clocks
sudo nvpmodel -m 2
sudo jetson_clocks

# Verify GPU is at max
GPU_FREQ=$(cat /sys/devices/platform/bus@0/17000000.gpu/devfreq/17000000.gpu/cur_freq)
echo "GPU frequency: $GPU_FREQ"
if [ "$GPU_FREQ" != "1020000000" ]; then
    echo "WARNING: GPU not locked to 1020 MHz, aborting"
    exit 1
fi

# Drop page cache for cold run
sync && echo 3 | sudo tee /proc/sys/vm/drop_caches > /dev/null

#DATA_DIR=${REALPATH}/test/data/positive/AS2603/realtime/
#DATA_DIR=${REALPATH}/test/data/QGXXOX250638/
#DATA_DIR=/home/william/realfshd/test/data/bench_QGXXOX250638
# DATA_DIR=/home/william/realfshd/test/data/sample_name
DATA_DIR=/home/william/realfshd/test/data/realdata
rm -f ${DATA_DIR}/*.log
rm -rf ${DATA_DIR}/*/*logs ${DATA_DIR}/*/slow5*
rm -rf ${DATA_DIR}/*/sam* ${DATA_DIR}/*/*logs ${DATA_DIR}/*/slow5*
"${REALPATH}/scripts/realfreq.sh" -n -m "${DATA_DIR}"