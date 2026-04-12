#!/bin/bash

set -x

# set -euo pipefail

# RAW INPUT
INPUT_UBAM=$1 #$/data/neynjo24/ubam_chunks/AS2603/chunk_001.bam /data/neynjo24/ubam/QGXHXX240243_fast_unaligned.bam #
OUTDIR=$2 #/data/neynjo24/results_realfshd/AS2603/chunk001
INPUT_BLOW5=$3 #/data/neynjo24/blow5/QGXHXX240243_reads_reduced.blow5 # $4
REBASECALLED_UBAM=$4

PREFIX="SAMPLE"

# Default Inputs
REALPATH=$(dirname "$(readlink -f "$0")")
REF=${REALPATH}/genome/hs1.fa
REGION_BED=${REALPATH}/inputs/d4z4_region.chm13.bed
MMI="${REALPATH}/inputs/hs1.mmi"

# Paths to tools
SAMTOOLS=samtools
MINIMAP2=minimap2
SLOW5TOOLS=slow5tools

# Parameters
MAPQ=30
MIN_READ_LENGTH=30000
MODE=map-ont

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# -------------------- Added logging setup --------------------
LOG_DIR=/home/william/realfshd/pipeline_logs
mkdir -p "$LOG_DIR"
PERF_LOG="$LOG_DIR/pre_step_perf.log"
TEGRASTATS_BIN=$(command -v tegrastats || true)
TIME_BIN=$(command -v /usr/bin/time || command -v time)

if [ ! -f "$PERF_LOG" ]; then
    echo -e "step\tstart_time\tend_time\twall_sec\tcpu_percent_time\ttegrastats_cpu_avg_percent\ttegrastats_cpu_max_percent\ttegrastats_gpu_avg_percent\ttegrastats_gpu_max_percent\tstatus" > "$PERF_LOG"
fi

run_logged() {
    local STEP_NAME="$1"
    shift
    local CMD_STR="$*"

    local START_EPOCH
    local END_EPOCH
    local START_HUMAN
    local END_HUMAN
    local WALL_SEC
    local TIME_TMP
    local TEGRA_TMP
    local TEGRA_PID=""
    local CPU_PCT="NA"
    local TG_CPU_AVG="NA"
    local TG_CPU_MAX="NA"
    local TG_GPU_AVG="NA"
    local TG_GPU_MAX="NA"
    local STATUS=0

    TIME_TMP=$(mktemp)
    TEGRA_TMP=$(mktemp)

    START_EPOCH=$(date +%s)
    START_HUMAN=$(date '+%Y-%m-%d %H:%M:%S')

    if [ -n "$TEGRASTATS_BIN" ]; then
        "$TEGRASTATS_BIN" --interval 1000 > "$TEGRA_TMP" 2>/dev/null &
        TEGRA_PID=$!
    fi

    /usr/bin/time -f 'CPU_PERCENT=%P' -o "$TIME_TMP" bash -c "$CMD_STR"
    STATUS=$?

    END_EPOCH=$(date +%s)
    END_HUMAN=$(date '+%Y-%m-%d %H:%M:%S')
    WALL_SEC=$((END_EPOCH - START_EPOCH))

    if [ -n "$TEGRA_PID" ]; then
        kill "$TEGRA_PID" 2>/dev/null || true
        wait "$TEGRA_PID" 2>/dev/null || true
    fi

    if grep -q '^CPU_PERCENT=' "$TIME_TMP" 2>/dev/null; then
        CPU_PCT=$(sed -n 's/^CPU_PERCENT=//p' "$TIME_TMP" | tail -n1)
    fi

    if [ -s "$TEGRA_TMP" ]; then
        read TG_CPU_AVG TG_CPU_MAX TG_GPU_AVG TG_GPU_MAX < <(
            gawk '
            {
                # -------- GPU: GR3D_FREQ 45% or GR3D_FREQ 0%@...
                if (match($0, /GR3D_FREQ[[:space:]]+([0-9]+)%/, g)) {
                    gpu = g[1] + 0
                    gpu_sum += gpu
                    gpu_n++
                    if (gpu > gpu_max) gpu_max = gpu
                }

                # -------- CPU: CPU [1%@1728,5%@1728,off,12%@1728]
                if (match($0, /CPU[[:space:]]+\[([^\]]+)\]/, c)) {
                    n = split(c[1], arr, ",")
                    local_sum = 0
                    local_n = 0
                    for (i = 1; i <= n; i++) {
                        if (match(arr[i], /([0-9]+)%/, p)) {
                            val = p[1] + 0
                            local_sum += val
                            local_n++
                            if (val > cpu_core_max) cpu_core_max = val
                        }
                    }
                    if (local_n > 0) {
                        sample_avg = local_sum / local_n
                        cpu_sample_sum += sample_avg
                        cpu_sample_n++
                        if (sample_avg > cpu_sample_max) cpu_sample_max = sample_avg
                    }
                }
            }
            END {
                cpu_avg = (cpu_sample_n > 0) ? cpu_sample_sum / cpu_sample_n : -1
                cpu_max = (cpu_sample_n > 0) ? cpu_sample_max : -1
                gpu_avg = (gpu_n > 0) ? gpu_sum / gpu_n : -1
                gpu_max_out = (gpu_n > 0) ? gpu_max : -1

                if (cpu_avg < 0) printf "NA "
                else printf "%.2f ", cpu_avg

                if (cpu_max < 0) printf "NA "
                else printf "%.2f ", cpu_max

                if (gpu_avg < 0) printf "NA "
                else printf "%.2f ", gpu_avg

                if (gpu_max_out < 0) printf "NA\n"
                else printf "%.2f\n", gpu_max_out
            }' "$TEGRA_TMP"
        )
    fi

    echo -e "${STEP_NAME}\t${START_HUMAN}\t${END_HUMAN}\t${WALL_SEC}\t${CPU_PCT}\t${TG_CPU_AVG}\t${TG_CPU_MAX}\t${TG_GPU_AVG}\t${TG_GPU_MAX}\t${STATUS}" >> "$PERF_LOG"

    rm -f "$TIME_TMP" "$TEGRA_TMP"

    return $STATUS
}
# ------------------ End added logging setup ------------------

# Convert to FASTQ
UBAM_BASE=$(basename $INPUT_UBAM)
INPUT_FASTQ="${OUTDIR}/${UBAM_BASE%.bam}.fastq"

run_logged "samtools_fastq" \
    "$SAMTOOLS fastq -TMM,ML \"$INPUT_UBAM\" > \"$INPUT_FASTQ\""

echo "minimap"
CHM13_BAM="${OUTDIR}/${PREFIX}_fast_chm13.bam"
SPLIT_PREFIX="${OUTDIR}/${PREFIX}_hs1_split_tmp"

echo "RUNNING: $MINIMAP2 -K 50M -x $MODE -a $MMI $INPUT_FASTQ --split-prefix $SPLIT_PREFIX -t 6 | $SAMTOOLS sort -o $CHM13_BAM" >> "$PERF_LOG"
run_logged "minimap2_sort" \
    "$MINIMAP2 -K 50M -x \"$MODE\" -a \"$MMI\" \"$INPUT_FASTQ\" --split-prefix \"$SPLIT_PREFIX\" -t 6 | $SAMTOOLS sort -o \"$CHM13_BAM\""

run_logged "samtools_index_chm13" \
    "$SAMTOOLS index \"$CHM13_BAM\""

echo "extracting reads mapping to the region of interest"
# Get d4z4 reads
D4Z4_ROI_BAM="${OUTDIR}/${PREFIX}_d4z4_reads_of_interest_fast.bam"
run_logged "extract_roi_sort" \
    "$SAMTOOLS view -L \"$REGION_BED\" -Sb \"$CHM13_BAM\" | $SAMTOOLS sort -o \"$D4Z4_ROI_BAM\""

run_logged "samtools_index_roi" \
    "$SAMTOOLS index \"$D4Z4_ROI_BAM\""

# -----------------------------------------------------------------------------------
#  Extract D4Z4 reads for rebasecalling with methylation
# -----------------------------------------------------------------------------------
d4z4_rid="${OUTDIR}/${PREFIX}_d4z4_rid.txt"

run_logged "extract_rid_from_roi" \
    "$SAMTOOLS view \"$D4Z4_ROI_BAM\" | cut -f1 | sort -u > \"$d4z4_rid\""

# Extract reads from uBAM
D4Z4_UBAM="${OUTDIR}/${PREFIX}_d4z4_reads_in_ubam.bam"
echo "extract d4z4 reads from ubam"
run_logged "extract_d4z4_ubam" \
    "$SAMTOOLS view -N \"$d4z4_rid\" -b \"$INPUT_UBAM\" > \"$D4Z4_UBAM\""

run_logged "samtools_index_d4z4_ubam" \
    "$SAMTOOLS index \"$D4Z4_UBAM\""

# Extract read IDs with and without pi tags from uBAM
# d4z4_rid="${OUTDIR}/${PREFIX}_d4z4_rid.txt}"

# Extract reads without pi tags
run_logged "extract_rid_without_pi" \
    "grep -vP 'pi:Z:\\K\\S+' <( $SAMTOOLS view \"$D4Z4_UBAM\" ) | cut -f1 | sort -u > \"$d4z4_rid\""

# Extract read with pi tags
run_logged "extract_rid_with_pi" \
    "grep -oP 'pi:Z:\\K\\S+' <( $SAMTOOLS view \"$D4Z4_UBAM\" ) | sort -u >> \"$d4z4_rid\""

run_logged "dedup_final_rid" \
    "sort -u \"$d4z4_rid\" > \"${OUTDIR}/d4z4_rid_final.txt\" && mv \"${OUTDIR}/d4z4_rid_final.txt\" \"$d4z4_rid\""

# Index INPUT BLOW5
run_logged "slow5_index_input" \
    "${SLOW5TOOLS} index \"${INPUT_BLOW5}\""

# After filtering reads with on D4Z4
echo "extract d4z4 reads from blow5"
EXTRACTED_BLOW5="${OUTDIR}/${PREFIX}_d4z4_reads.blow5"
run_logged "slow5_get_d4z4" \
    "${SLOW5TOOLS} get \"${INPUT_BLOW5}\" --list \"${d4z4_rid}\" -o \"${EXTRACTED_BLOW5}\""

run_logged "slow5_index_extracted" \
    "${SLOW5TOOLS} index \"${EXTRACTED_BLOW5}\""

# rebasecall with SUP
echo "rebasecall with SUP"
run_logged "slow5_dorado_sup" \
    "${SLOW5_DORADO} basecaller -x cuda:all --models-directory ${DORADO_MODEL_PATH} ${DORADO_MODEL_PATH}/${DORADO_MODEL_SUP} \"${EXTRACTED_BLOW5}\" --modified-bases 5mCG_5hmCG > \"${REBASECALLED_UBAM}\""

# Remove blow5
run_logged "cleanup" \
    "rm -f \"${EXTRACTED_BLOW5}\" \"${EXTRACTED_BLOW5}.idx\" \"${INPUT_FASTQ}\" \"${D4Z4_ROI_BAM}\" \"${D4Z4_ROI_BAM}.bai\""
# # rm "$d4z4_rid"