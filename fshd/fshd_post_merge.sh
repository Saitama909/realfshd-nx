#!/bin/bash

# ----------------------------------------------------------------
# Realtime stop analysis
# 1. Merge myop panel BAM
# 2. merge feature summary
# 3. run fshd1 status generation
# 4. run realtime_stop.sh

OUTDIR=$1
PAR_DIR=$(dirname "$OUTDIR")
SAMPLE_NAME=SAMPLE
MERGED_MYOP_BAM="${PAR_DIR}/merged/${SAMPLE_NAME}_myop_panel.bam"
CHUNK_MYOP_BAM="${OUTDIR}/${SAMPLE_NAME}_myop_panel.bam"

MERGED_FEATURES_SUMMARY_TSV="${PAR_DIR}/merged/${SAMPLE_NAME}_mapped_features_summary.tsv"
CHUNK_FEATURES_SUMMARY_TSV="${OUTDIR}/${SAMPLE_NAME}_mapped_features_summary.tsv"

REALPATH=$(dirname "$(readlink -f "$0")")

echo "Merging myopathy panel BAM to: ${MERGED_MYOP_BAM}"
echo "Merging features summary to: ${MERGED_FEATURES_SUMMARY_TSV}"

# -------------------- Added logging setup --------------------
LOG_DIR="/home/william/realfshd/pipeline_logs"
mkdir -p "$LOG_DIR"
PERF_LOG="${LOG_DIR}/post_step_perf.log"

TEGRASTATS_BIN=$(command -v tegrastats || true)
GAWK_BIN=$(command -v gawk || true)

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

    /usr/bin/time -f 'CPU_PERCENT=%P' -o "$TIME_TMP" bash -o pipefail -c "$CMD_STR"
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

    if [ -n "$GAWK_BIN" ] && [ -s "$TEGRA_TMP" ]; then
        read TG_CPU_AVG TG_CPU_MAX TG_GPU_AVG TG_GPU_MAX < <(
            gawk '
            {
                if (match($0, /GR3D_FREQ[[:space:]]+([0-9]+)%/, g)) {
                    gpu = g[1] + 0
                    gpu_sum += gpu
                    gpu_n++
                    if (gpu > gpu_max) gpu_max = gpu
                }

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

    echo -e "${STEP_NAME}\t${START_HUMAN}\t${END_HUMAN}\t${WALL_SEC}\t${CPU_PCT}\t${TG_CPU_AVG}\t${TG_CPU_MAX}\t${TG_GPU_AVG}\t${TG_GPU_MAX}\t${STATUS}" | tee -a "$PERF_LOG"

    rm -f "$TIME_TMP" "$TEGRA_TMP"
    return $STATUS
}
# ------------------ End added logging setup ------------------

run_logged "prepare_merged_dir" \
    "mkdir -p \"${PAR_DIR}/merged\""

if [ ! -f "${MERGED_MYOP_BAM}" ] || [ ! -f "${MERGED_FEATURES_SUMMARY_TSV}" ]; then
    run_logged "init_copy_myop_bam" \
        "cp \"${CHUNK_MYOP_BAM}\" \"${MERGED_MYOP_BAM}\""

    run_logged "init_index_myop_bam" \
        "samtools index \"${MERGED_MYOP_BAM}\""

    run_logged "init_copy_features_summary" \
        "cp \"${CHUNK_FEATURES_SUMMARY_TSV}\" \"${MERGED_FEATURES_SUMMARY_TSV}\""
else
    run_logged "append_features_summary" \
        "tail -n +2 \"${CHUNK_FEATURES_SUMMARY_TSV}\" >> \"${MERGED_FEATURES_SUMMARY_TSV}\""

    tmp_bam="${MERGED_MYOP_BAM%.bam}.tmp.bam"
    run_logged "merge_myop_bams" \
        "samtools merge -f \"$tmp_bam\" \"$MERGED_MYOP_BAM\" \"$CHUNK_MYOP_BAM\" && mv \"$tmp_bam\" \"$MERGED_MYOP_BAM\""

    run_logged "index_merged_myop_bam" \
        "samtools index \"${MERGED_MYOP_BAM}\""
fi

run_logged "extract_fshd1_status" \
    "source \"${REALPATH}/venv3/bin/activate\" && PYTHONPATH=\"${REALPATH}/helper\" python3 -c '
from methylation_summary import extract_fshd1_status
import pandas as pd
summary_df = pd.read_csv(\"${MERGED_FEATURES_SUMMARY_TSV}\", sep=\"\t\")
extract_fshd1_status(summary_df, \"${PAR_DIR}/merged/${SAMPLE_NAME}_fshd1_status_counts.tsv\", \"${REALPATH}/inputs/fshd1_status_template.tsv\")
'"

run_logged "realtime_stop" \
    "\"${REALPATH}/realtime_stop.sh\" \
    \"${PAR_DIR}/merged/${SAMPLE_NAME}_fshd1_status_counts.tsv\" \
    \"${PAR_DIR}/merged/${SAMPLE_NAME}_mapped_features_summary.tsv\" \
    \"${PAR_DIR}/merged/${SAMPLE_NAME}_myop_panel.bam\""

OUTDIR_COUNTER="${OUTDIR}/merged"
run_logged "copy_merged_record" \
    "mkdir -p \"${OUTDIR_COUNTER}\" && cp \"${PAR_DIR}/merged\"/* \"${OUTDIR_COUNTER}/.\""