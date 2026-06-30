#!/bin/bash

# Default values
TIME_LOG="realfreq_pipeline_time.log"
TMP_FAILED="realfreq_failed_list.log"
MAX_PROC=1
CHECK_TOOLS=false

# Colour codes
YELLOW="\e[33m"
RED="\e[31m"
NORMAL="\033[0;39m"

# terminate script
die() { echo -e ${RED}"$1"${NORMAL} >&2; echo; exit 1; }
failed() { echo -e ${RED}"$1"${NORMAL} >&2; echo; echo "$2" >> "$TMP_FAILED"; exit 1; }

## Handle flags
while getopts "l:f:p:c" o; do
    case "${o}" in
        l) TIME_LOG=${OPTARG} ;;
        f) TMP_FAILED=${OPTARG} ;;
        p) MAX_PROC=${OPTARG} ;;
        c) CHECK_TOOLS=true ;;
        *) echo "[pipeline.sh] Incorrect args"; exit 1 ;;
    esac
done
shift $((OPTIND-1))

#====================User Start - Tool check====================

[ -z "${SLOW5TOOLS}" ]            && export SLOW5TOOLS=slow5tools
[ -z "${BLUECRAB}" ]              && export BLUECRAB=blue-crab
[ -z "${SLOW5_DORADO}" ]          && export SLOW5_DORADO=slow5-dorado
[ -z "${FSHD_PRE}" ]              && export FSHD_PRE=pre_mapping_n_rebasecall.sh
[ -z "${FSHD_PIPELINE_WRAPPER}" ] && export FSHD_PIPELINE_WRAPPER=fshd_pipeline_wrapper.sh
[ -z "${FSHD_POST}" ]             && export FSHD_POST=fshd_post_merge.sh
[ -z "${FSHDTOOLS}" ]             && export FSHDTOOLS=fshdtools

"$SLOW5TOOLS" --version &> /dev/null   || die "[pipeline.sh] slow5tools not found. Add to PATH or export SLOW5TOOLS=/path/to/slow5tools. Exiting."
"$BLUECRAB" --version &> /dev/null     || die "[pipeline.sh] bluecrab not found. Add to PATH or export BLUECRAB=/path/to/bluecrab. Exiting."
"$SLOW5_DORADO" --version &> /dev/null || die "[pipeline.sh] slow5-dorado not found. Add to PATH or export SLOW5_DORADO=/path/to/slow5-dorado. Exiting."
"$FSHDTOOLS" --version &> /dev/null    || die "[pipeline.sh] fshdtools not found. Add to PATH or export FSHDTOOLS=/path/to/fshdtools. Exiting."

[ -z "${DORADO_MODEL_PATH}" ] && die "[pipeline.sh] MODEL_PATH not set. export MODEL_PATH=/path/to/model. Exiting."
[ -z "${DORADO_MODEL_FAST}" ] && die "[pipeline.sh] DORADO_MODEL_FAST not set. export DORADO_MODEL=dorado_model. Exiting."
[ -z "${DORADO_MODEL_SUP}" ]  && die "[pipeline.sh] DORADO_MODEL_SUP not set. export DORADO_MODEL=dorado_model. Exiting."

#====================User End - Tool check====================

if [ "$CHECK_TOOLS" = true ]; then
    echo "[pipeline.sh] Tool check successful"
    exit 0
fi

echo "[pipeline.sh] Starting pipeline with $MAX_PROC max processes"
counter=0

while read -r PIPELINE_INPUT
do
(
    #=======================User Start - Pipeline=======================

    P5_FILEPATH="$PIPELINE_INPUT"
    P5_DIR=${P5_FILEPATH%/*}
    PARENT_DIR=${P5_DIR%/*}
    P5_FILENAME=$(basename "$P5_FILEPATH")
    P5_PREFIX=${P5_FILENAME%.*}

    if [[ "$P5_DIR" =~ .*pod5_pass.* ]]; then
        SLOW5_DIR=$(echo "$P5_DIR" | sed 's/pod5_pass/slow5_pass/g')
        SAM_DIR=$(echo "$P5_DIR"   | sed 's/pod5_pass/sam_pass/g')
        LOG_DIR=$(echo "$P5_DIR"   | sed 's/pod5_pass/pass_logs/g')
    elif [[ "$P5_DIR" =~ .*pod5_fail.* ]]; then
        SLOW5_DIR=$(echo "$P5_DIR" | sed 's/pod5_fail/slow5_fail/g')
        SAM_DIR=$(echo "$P5_DIR"   | sed 's/pod5_fail/sam_fail/g')
        LOG_DIR=$(echo "$P5_DIR"   | sed 's/pod5_fail/fail_logs/g')
    elif [[ "$P5_DIR" =~ .*pod5_skip.* ]]; then
        SLOW5_DIR=$(echo "$P5_DIR" | sed 's/pod5_skip/slow5_skip/g')
        SAM_DIR=$(echo "$P5_DIR"   | sed 's/pod5_skip/sam_skip/g')
        LOG_DIR=$(echo "$P5_DIR"   | sed 's/pod5_skip/skip_logs/g')
    else
        SLOW5_DIR="$PARENT_DIR/slow5"
        SAM_DIR="$PARENT_DIR/sam"
        LOG_DIR="$PARENT_DIR/logs"
    fi

    if [ -z "$SLOW5_DIR" ] || [ -z "$LOG_DIR" ]; then
        SLOW5_DIR="$PARENT_DIR/slow5"
        SAM_DIR="$PARENT_DIR/sam"
        LOG_DIR="$PARENT_DIR/logs"
    fi

    test -d "$SLOW5_DIR" || { mkdir -p "$SLOW5_DIR"; echo "[pipeline.sh] Created $SLOW5_DIR/. Converted SLOW5 files will be here."; }
    test -d "$LOG_DIR"   || { mkdir -p "$LOG_DIR";   echo "[pipeline.sh] Created $LOG_DIR/. individual logs for each conversion will be here."; }
    test -d "$SAM_DIR"   || { mkdir -p "$SAM_DIR";   echo "[pipeline.sh] Created $SAM_DIR/. SAM files will be here."; }

    SLOW5_FILEPATH="$SLOW5_DIR/$P5_PREFIX.blow5"
    SLOW5_LONG_FILEPATH="$SLOW5_DIR/$P5_PREFIX.long.blow5"
    UNALN_SAM_FILEPATH="$SAM_DIR/$P5_PREFIX.fast.unaln.sam"
    FASTQ_FILEPATH="$SAM_DIR/$P5_PREFIX.fastq"
    FSHD_OUTPATH="$SAM_DIR/$P5_PREFIX.fshd.out"
    REBASECALLED_UBAM="$SAM_DIR/$P5_PREFIX.sup.unsorted.bam"
    LOG_FILEPATH="$LOG_DIR/$P5_PREFIX.log"

    TEGRASTATS_BIN=$(command -v tegrastats || true)
    GAWK_BIN=$(command -v gawk || true)

    run_step() {
        local STEP_NAME="$1"
        local CMD_STR="$2"
        local t_start t_end status

        local TIME_TMP TEGRA_TMP TEGRA_PID=""
        local CPU_PCT="NA"
        local RAM_AVG_MB="NA"
        local CPU_FREQ_MIN="NA"   # ADDED: minimum CPU frequency during step (MHz) — drops if throttling

        TIME_TMP=$(mktemp)
        TEGRA_TMP=$(mktemp)

        t_start=$(date '+%Y-%m-%d %H:%M:%S')
        echo "[pipeline.sh] START ${STEP_NAME} ${P5_FILEPATH} at ${t_start}" | tee -a "$LOG_FILEPATH"

        if [ -n "$TEGRASTATS_BIN" ]; then
            "$TEGRASTATS_BIN" --interval 1000 > "$TEGRA_TMP" 2>/dev/null &
            TEGRA_PID=$!
        fi

        set +e
        /usr/bin/time -f 'CPU_PERCENT=%P' -o "$TIME_TMP" bash -o pipefail -c "$CMD_STR" >> "$LOG_FILEPATH" 2>&1
        status=$?
        set -e

        t_end=$(date '+%Y-%m-%d %H:%M:%S')

        if [ -n "$TEGRA_PID" ]; then
            kill "$TEGRA_PID" 2>/dev/null || true
            wait "$TEGRA_PID" 2>/dev/null || true
        fi

        # Extract CPU% from /usr/bin/time output
        if grep -q '^CPU_PERCENT=' "$TIME_TMP" 2>/dev/null; then
            CPU_PCT=$(sed -n 's/^CPU_PERCENT=//p' "$TIME_TMP" | tail -n1)
        fi

        # Extract average RAM used (MB) from tegrastats output
        # tegrastats reports RAM like: RAM 3012/7772MB
        if [ -n "$GAWK_BIN" ] && [ -s "$TEGRA_TMP" ]; then
            RAM_AVG_MB=$(gawk '
            {
                if (match($0, /RAM[[:space:]]+([0-9]+)\/([0-9]+)MB/, r)) {
                    sum += r[1]+0
                    n++
                }
            }
            END {
                if (n > 0) printf "%.0f", sum/n
                else print "NA"
            }' "$TEGRA_TMP")

            # ADDED: extract minimum CPU frequency during step
            # tegrastats reports CPU like: CPU [45%@1728,50%@1728,off,12%@1728]
            # frequency drops below 1728 MHz indicates thermal throttling
            CPU_FREQ_MIN=$(gawk '
            {
                if (match($0, /CPU[[:space:]]+\[([^\]]+)\]/, c)) {
                    n = split(c[1], arr, ",")
                    for (i = 1; i <= n; i++) {
                        if (match(arr[i], /@([0-9]+)/, f)) {
                            freq = f[1] + 0
                            if (min == 0 || freq < min) min = freq
                        }
                    }
                }
            }
            END { print (min > 0 ? min : "NA") }' "$TEGRA_TMP")
        fi

        rm -f "$TIME_TMP" "$TEGRA_TMP"

        if [ "$status" -eq 0 ]; then
            echo "[pipeline.sh] SUCCESS ${STEP_NAME} ${P5_FILEPATH} at ${t_end}" | tee -a "$LOG_FILEPATH"
            echo -e "${P5_FILEPATH}\t${STEP_NAME}\t${t_start}\t${t_end}\tSUCCESS\t${status}\t${CPU_PCT}\t${RAM_AVG_MB}\t${CPU_FREQ_MIN}" >> "${TIME_LOG}"
        else
            echo "[pipeline.sh] FAIL ${STEP_NAME} ${P5_FILEPATH} at ${t_end} (exit=${status})" | tee -a "$LOG_FILEPATH"
            echo -e "${P5_FILEPATH}\t${STEP_NAME}\t${t_start}\t${t_end}\tFAIL\t${status}\t${CPU_PCT}\t${RAM_AVG_MB}\t${CPU_FREQ_MIN}" >> "${TIME_LOG}"
            failed "${STEP_NAME} failed. check $LOG_FILEPATH" "$P5_FILEPATH"
        fi
    }
    run_step "p2s" \
        "\"${BLUECRAB}\" p2s -p1 -c zstd -s svb-zd \"$P5_FILEPATH\" -o \"$SLOW5_FILEPATH\""

    # run_step "p2s" \
    #     "\"${BLUECRAB}\" p2s -p1 \"$P5_FILEPATH\" -o \"$SLOW5_FILEPATH\""

    run_step "long" \
        "\"${FSHDTOOLS}\" long \"$SLOW5_FILEPATH\" -o \"$SLOW5_LONG_FILEPATH\""

    run_step "fastbasecall" \
        "${SLOW5_DORADO} basecaller -x cuda:0 --slow5_batchsize 64 --models-directory ${DORADO_MODEL_PATH} ${DORADO_MODEL_PATH}/${DORADO_MODEL_FAST} ${SLOW5_LONG_FILEPATH} > ${UNALN_SAM_FILEPATH}"

    run_step "pre" \
        "\"${FSHD_PRE}\" \"$UNALN_SAM_FILEPATH\" \"$SAM_DIR\" \"$SLOW5_LONG_FILEPATH\" \"$REBASECALLED_UBAM\""

    run_step "fshd" \
        "\"${FSHD_PIPELINE_WRAPPER}\" \"$REBASECALLED_UBAM\" \"${FSHD_OUTPATH}/\""

    run_step "post" \
        "\"${FSHD_POST}\" \"${FSHD_OUTPATH}/\""

    echo "[pipeline.sh] Done $P5_FILEPATH" | tee -a "$LOG_FILEPATH"

    PIPELINE_OUTPUT="$REBASECALLED_UBAM"
    #=======================User End - Pipeline=======================

    echo "realfreq-pipeline-output:$PIPELINE_OUTPUT"
)&
    ((counter++))
    if [ "$counter" -ge "$MAX_PROC" ]; then
        echo "[pipeline.sh] Waiting for $counter jobs to finish."
        wait
        counter=0
    fi
done

wait