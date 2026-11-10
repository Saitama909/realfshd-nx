#!/bin/bash

# module unload samtools
# module unload bedtools
# module load samtools/1.22
# module load bedtools/2.31.0

STATUS_TSV=$1
SUMMARY_TSV=$2 # "${OUTDIR}/${PREFIX}_mapped_features_summary.tsv"
MYOP_BAM=$3 # "${OUTDIR}/${PREFIX}_myop_panel.bam"

REALPATH=$(dirname "$(readlink -f "$0")")
MYOP_PANEL_BED=${REALPATH}/inputs/myopathy_panel_draft_2.gene_targets.chm13.bed
SAMPLE_DIR=$(dirname "$STATUS_TSV")
# MIN_COMPLETE4QA=5

# Check if summary tsv exist
if [[ ! -f "$STATUS_TSV" ]]; then
  echo "Error: file $STATUS_TSV not found"
  exit 1
fi

# Check if summary tsv exist
if [[ ! -f "$SUMMARY_TSV" ]]; then
  echo "Error: file $SUMMARY_TSV not found"
  exit 1
fi

# Check if myop bam exist
if [[ ! -f "$MYOP_BAM" ]]; then
  echo "Error: file $MYOP_BAM not found"
  exit 1
fi

# ---------------------------------------------------------------------------------------------
# Logging setup (ADDED)
# ---------------------------------------------------------------------------------------------
LOG_DIR="/home/william/realfshd/pipeline_logs"
mkdir -p "$LOG_DIR"
PERF_LOG="${LOG_DIR}/realtime_stop_step_perf.log"

TEGRASTATS_BIN=$(command -v tegrastats || true)
GAWK_BIN=$(command -v gawk || true)

if [ ! -f "$PERF_LOG" ]; then
    echo -e "step\tstart_time\tend_time\twall_sec\tcpu_percent\ttegrastats_cpu_avg\ttegrastats_cpu_max\ttegrastats_gpu_avg\ttegrastats_gpu_max\tstatus" > "$PERF_LOG"
fi

run_logged() {
    local STEP_NAME="$1"
    shift
    local CMD_STR="$*"

    local TIME_TMP TEGRA_TMP TEGRA_PID=""
    local CPU_PCT="NA"
    local TG_CPU_AVG="NA" TG_CPU_MAX="NA"
    local TG_GPU_AVG="NA" TG_GPU_MAX="NA"
    local STATUS=0

    TIME_TMP=$(mktemp)
    TEGRA_TMP=$(mktemp)

    local START_EPOCH START_HUMAN END_EPOCH END_HUMAN WALL_SEC
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
                    gpu = g[1] + 0; gpu_sum += gpu; gpu_n++
                    if (gpu > gpu_max) gpu_max = gpu
                }
                if (match($0, /CPU[[:space:]]+\[([^\]]+)\]/, c)) {
                    n = split(c[1], arr, ",")
                    local_sum = 0; local_n = 0
                    for (i = 1; i <= n; i++) {
                        if (match(arr[i], /([0-9]+)%/, p)) {
                            val = p[1] + 0; local_sum += val; local_n++
                            if (val > cpu_core_max) cpu_core_max = val
                        }
                    }
                    if (local_n > 0) {
                        sample_avg = local_sum / local_n
                        cpu_sample_sum += sample_avg; cpu_sample_n++
                        if (sample_avg > cpu_sample_max) cpu_sample_max = sample_avg
                    }
                }
            }
            END {
                cpu_avg = (cpu_sample_n > 0) ? cpu_sample_sum / cpu_sample_n : -1
                cpu_max = (cpu_sample_n > 0) ? cpu_sample_max : -1
                gpu_avg = (gpu_n > 0) ? gpu_sum / gpu_n : -1
                gpu_max_out = (gpu_n > 0) ? gpu_max : -1
                printf "%s %s %s %s\n",
                    (cpu_avg < 0 ? "NA" : sprintf("%.2f", cpu_avg)),
                    (cpu_max < 0 ? "NA" : sprintf("%.2f", cpu_max)),
                    (gpu_avg < 0 ? "NA" : sprintf("%.2f", gpu_avg)),
                    (gpu_max_out < 0 ? "NA" : sprintf("%.2f", gpu_max_out))
            }' "$TEGRA_TMP"
        )
    fi

    local LOG_LINE
    LOG_LINE=$(printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
        "$STEP_NAME" "$START_HUMAN" "$END_HUMAN" "$WALL_SEC" "$CPU_PCT" \
        "$TG_CPU_AVG" "$TG_CPU_MAX" "$TG_GPU_AVG" "$TG_GPU_MAX" "$STATUS")
    printf '%s\n' "$LOG_LINE" >> "$PERF_LOG"
    printf '%s\n' "$LOG_LINE" >&2

    rm -f "$TIME_TMP" "$TEGRA_TMP"
    return $STATUS
}
# ---------------------------------------------------------------------------------------------
# End logging setup
# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Set thresholds
# ---------------------------------------------------------------------------------------------
POSITIVE_THRESH=3 # Positive FSHD
NEGATIVE_THRESH=5 # Negative FSHD
COVERAGE_THRESH=30 # Threshold for gene to stop being in read fish from panel
COVERAGE_THRESH_CHRX=30
LOD_THRESH=0.95 # Limit of detection threshold for FSHD

# ---------------------------------------------------------------------------------------------
# Positive case - Stop sequencing
# ---------------------------------------------------------------------------------------------

positive_status=false
pos_out=$(run_logged "positive_fshd_check" \
    "awk -v thresh='$POSITIVE_THRESH' -F'\t' '\$1==\"Diagnostic/Positive\" && \$7>=thresh { print \"true\" }' '$STATUS_TSV'")
[ -n "$pos_out" ] && read -r positive_status <<< "$pos_out"

# Check awk exit code
# Stop the script if it's enough
if [ "$positive_status" = "true" ]; then
    readlabel_col=$(awk -F'\t' 'NR==1 { for (i=1; i<=NF; i++) if ($i=="ReadLabel") print i }' $SUMMARY_TSV)
	copies_col=$(awk -F'\t' 'NR==1 { for (i=1; i<=NF; i++) if ($i=="MappedEstimatedCopies") print i }' $SUMMARY_TSV)

	# Check if all the complete reads have the same floored length
	# And if not, check how many reads each length has, and if there are more than 2 different lengths, stop the script and ask for manual review
	echo "Checking if all complete reads have the same length..."
	num_copies=$(awk -F'\t' -v rl_col="$readlabel_col" -v cp_col="$copies_col" 'NR>1 && $rl_col=="Complete 4qA" { print int($cp_col) }' $SUMMARY_TSV | sort | uniq -c | sed 's/^[ \t]*//')

	# IF there is more than 1 different length, check how many reads in each length, if there are more than >=3 reads each length, stop the script
	num_lengths=$(echo "$num_copies" | wc -l)
	if [ "$num_lengths" -gt 1 ]; then
		num_reads=$(echo "$num_copies" | sort | tail -n1 | cut -f1 -d" ")
		if [ "$num_reads" -ge 3 ]; then
			echo "✅ Positive FSHD1 case, STOP SEQUENCING"
			exit 0
		else
			echo "❌ Warning: More than 1 different contracted length found, but not enough reads in each length, continue sequencing"
		fi
	else
		echo "✅ Positive FSHD1 case, STOP SEQUENCING"
		exit 0
	fi
fi

# ---------------------------------------------------------------------------------------------
# Easy negative case (1x 4qA with more than 10 copies and 1x 4qB)
# ---------------------------------------------------------------------------------------------
neg_out=$(run_logged "negative_fshd_check" \
    "awk -F'\t' '\$1 == \"Diagnostic/Negative\" { count++; if (\$NF >= 3) ok++ } END { if (count == 2 && ok == 2) print \"True\", \"✅ Have a non-permissive 4qA and 4qB reads with ≥3 reads each\"; else print \"False\", \"❌ One or both Diagnostic/Negative rows have <3 reads\" }' '$STATUS_TSV'")
read negative_status message <<< "$neg_out"

echo "$message"

# ---------------------------------------------------------------------------------------------
# Easy negative case 2 *only run if there's a complete reads for 4q*
# either:
# 2x 4qA with more than 10 copies on MappedEstimatedCopies and both >= NEGATIVE_THRESH
# 2x 4qB with different MappedEstimatedCopies and both >= NEGATIVE_THRESH
# Check two different lengths by getting MappedEstimatedCopies for 4qA/4qB row from SUMMARY_TSV
# Floor the MappedEstimatedCopies to integer for counting purpose
# If both lengths have >= NEGATIVE_THRESH, then stop sequencing
# ---------------------------------------------------------------------------------------------
run_logged "negative_fshd_check_2" \
    "grep -q 'Complete 4q' '$SUMMARY_TSV' && echo 'complete 4q found' || echo 'no complete 4q'"

grep -q "Complete 4q" $SUMMARY_TSV && awk -F'\t' -v threshold="$NEGATIVE_THRESH" '
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i == "ReadLabel") rcol = i
        if ($i == "MappedEstimatedCopies") ccol = i
    }
    next
}
{
    label = $rcol
    copies = int($ccol)

    if (label == "Complete 4qA") {
        if (copies > 11)
            count[label,copies]++
    }
    else if (label == "Complete 4qB") {
        count[label,copies]++
    }
}

END {
    # Count distinct values
    for (key in count) {
        split(key, parts, SUBSEP)
        lab = parts[1]
        cop = parts[2]
        distinct[lab]++
        summary[lab] = summary[lab] " " cop ":" count[key]
    }

    ### CHECK 4qA ###
    okA = 1
    if (distinct["Complete 4qA"] < 2)
        okA = 0
    else {
        for (key in count) {
            split(key, parts, SUBSEP)
            if (parts[1] == "Complete 4qA" && count[key] < threshold)
                okA = 0
        }
    }

    ### CHECK 4qB ###
    okB = 1
    if (distinct["Complete 4qB"] < 2)
        okB = 0
    else {
        for (key in count) {
            split(key, parts, SUBSEP)
            if (parts[1] == "Complete 4qB" && count[key] < threshold)
                okB = 0
        }
    }

    print "4qA distinct:", distinct["Complete 4qA"], "→", summary["Complete 4qA"]
    print "4qB distinct:", distinct["Complete 4qB"], "→", summary["Complete 4qB"]

    print (okA ? "Complete 4qA OK" : "Complete 4qA FAIL")
    print (okB ? "Complete 4qB OK" : "Complete 4qB FAIL")
}' "$SUMMARY_TSV"

# ---------------------------------------------------------------------------------------------
# Coverage check for myopathy panel
# ---------------------------------------------------------------------------------------------
# Only check coverage if not all genes have complete coverage
# function triggered when all autosomal chromsome reach coverage threshold but some genes on chrX haven't reach the chrX coverage threshold
# Check if chrX coverage is on average half autosomal chr, chrX threshold become 15x, otherwise keep 30x
check_sex_coverage() {
    PANEL_COV_TXT="${SAMPLE_DIR}/all_myopathy_panel_chm13_coverage.txt"
    PANEL_COV_BED="${SAMPLE_DIR}/all_myopathy_panel_chm13_coverage.bed"

    run_logged "check_sex_coverage_depth" \
        "samtools depth -@ 6 -a -b '$MYOP_PANEL_BED' '$MYOP_BAM' | awk -F'\t' '{print \$1\"\t\"\$2\"\t\"\$2\"\t\"\$3}' > '$PANEL_COV_BED'"

    NOY_PANEL="${SAMPLE_DIR}/myop_panel_noY.bed"
    grep -v "^chrY" "$MYOP_PANEL_BED" > "$NOY_PANEL"

    run_logged "check_sex_coverage_intersect" \
        "bedtools intersect -a '$NOY_PANEL' -b '$PANEL_COV_BED' -wa -wb > ${SAMPLE_DIR}/all_panel_cov_intersect.bed"
    
    # Check average coverage of chrX and autosomal chr
    read sex coverage_thresh_chrX avgX avgAuto diff < <(
    awk -F'\t' '
    $1 == "chrX" { sumX += $8; countX++ }
    $1 != "chrX" { sumAuto += $8; countAuto++ }
    END {
        avgX = (countX > 0) ? sumX / countX : 0
        avgAuto = (countAuto > 0) ? sumAuto / countAuto : 0
        diff = (avgAuto / 2) - avgX
        if (diff < 0) diff = -diff
        if (diff <= 5) {
            print "male", 15, avgX, avgAuto, diff
        } else {
            print "female", 30, avgX, avgAuto, diff
        }
    }' ${SAMPLE_DIR}/all_panel_cov_intersect.bed
    )

    echo "$sex $coverage_thresh_chrX $avgX $avgAuto $diff"
}

TMP_PANEL="${SAMPLE_DIR}/myopathy_panel_chm13_temp.bed"
COVERAGE_OUTPUT="${SAMPLE_DIR}/myopathy_panel_chm13_coverage.txt"
output="${SAMPLE_DIR}/genes_complete_coverage.txt"
coverage_function() {
    echo "Checking coverage for myopathy panel genes..."
    # Duplicate myop bed file for reducing purpose once coverage reached
    # Copy panel if tmp doesn't exist
    if [ ! -f $TMP_PANEL ]; then
        cp $MYOP_PANEL_BED $TMP_PANEL
        > $COVERAGE_OUTPUT
        > ${SAMPLE_DIR}/genes_complete_coverage.txt
    fi

    echo "COVERAGE"
    # Samtools depth to get coverage over myopathy panel genes
    # Currently col 2 is coordinate and col 3 is depth
    # want to make it col 2 and col 3 is coordinate (duplicated col2)
    BED_COV="${SAMPLE_DIR}/myop_panel_chm13_coverage.bed"
    run_logged "samtools_depth_panel" \
        "samtools depth -@ 6 -a -b '$TMP_PANEL' '$MYOP_BAM' | awk -F'\t' '{print \$1\"\t\"\$2\"\t\"\$2\"\t\"\$3}' > '$BED_COV'"

    echo "INTERSECT"
    # get intersected regions of coverage and panel
    run_logged "bedtools_intersect_panel" \
        "bedtools intersect -a $TMP_PANEL -b $BED_COV -wa -wb > ${SAMPLE_DIR}/panel_cov_intersect.bed"

    if awk -F'\t' '
    NR > 1 && $1 != "chrX" && $1 != "chrY" { exit 1 }
    ' ${SAMPLE_DIR}/panel_cov_intersect.bed
    then
        echo "All autosomal genes have meet threshold, checking chrX genes"
        read sex COVERAGE_THRESH_CHRX avgX avgAuto diff < <(check_sex_coverage)
        echo "Sex: $sex, chrX coverage threshold: $COVERAGE_THRESH_CHRX, avgX: $avgX, avgAuto: $avgAuto, diff: $diff"
        # COVERAGE_THRESH_CHRX=$(check_sex_coverage)
    fi

    run_logged "awk_coverage_summary" \
        "awk -v cov_thresh='$COVERAGE_THRESH' \
        -v chrX_thresh='$COVERAGE_THRESH_CHRX' \
        -F'\t' '
    {
        chrom = \$1
        gene = \$4
        cov  = \$8

        sum_cov[gene] += cov
        base_count[gene]++
        gene_chr[gene] = chrom

        if (chrom == \"chrX\" && cov < chrX_thresh)
            under[gene]++
        else if (chrom != \"chrX\" && cov < cov_thresh)
            under[gene]++

    }

    END {
        printf \"%-20s\t%-15s\t%-10s\t%-15s\t%-18s\t%-15s\t%s\n\",
            \"Gene\", \"Avg_Coverage\", \"Length\", \"Total_Coverage\",
            \"Bases_Under_Threshold\", \"Pct_Bases_Above\", \"Pass_Threshold\"

        for (g in sum_cov) {
            total = base_count[g]
            under_count = (g in under ? under[g] : 0)
            above_pct = 1 - (under_count / total)
            avg = sum_cov[g] / total

            if (gene_chr[g] == \"chrX\")
                threshold = chrX_thresh
            else
                threshold = cov_thresh

            pass = (avg >= threshold && above_pct >= 0.95 ? \"YES\" : \"NO\")

            printf \"%-20s\t%-15.2f\t%-10d\t%-15d\t%-18d\t%-15.3f\t%s\n\",
                g, avg, total, sum_cov[g], under_count, above_pct, pass
        }
    }' '${SAMPLE_DIR}/panel_cov_intersect.bed' > '$output'"

    # Single-pass exact-match removal of genes that passed coverage
    PASSFILE="${SAMPLE_DIR}/genes_complete_coverage.txt"
    if awk -F'\t' -v passfile="$PASSFILE" '
        FILENAME == passfile {
            if (FNR == 1) next                        # skip printf header
            g = $1; gsub(/[[:space:]]+$/, "", g)      # strip printf padding
            p = $7; gsub(/[[:space:]]+$/, "", p)
            if (p == "YES") drop[g] = 1
            next
        }
        !($4 in drop)
    ' "$PASSFILE" "$TMP_PANEL" > "${TMP_PANEL}.tmp"
    then
        mv "${TMP_PANEL}.tmp" "$TMP_PANEL"
    else
        echo "⚠️  Gene-removal awk failed; leaving $TMP_PANEL unchanged"
        rm -f "${TMP_PANEL}.tmp"
    fi
}
    # if (print_header == 1) {
    #     }

# if lines in tmp panel = lines in MYOP PANEL
# then all genes have complete coverage
if [ -f ${TMP_PANEL} ]; then
    num_genes_to_complete=$(wc -l < <(grep -v "chrY" ${TMP_PANEL}))
else
    num_genes_to_complete=$(wc -l < <(grep -v "chrY" ${MYOP_PANEL_BED}))
fi
# if [ -f ${SAMPLE_DIR}/genes_complete_coverage.txt ]; then
#     num_genes_complete=$(wc -l < ${SAMPLE_DIR}/genes_complete_coverage.txt)
# else
#     num_genes_complete=0
# fi
# num_genes_complete=$(wc -l < ${SAMPLE_DIR}/genes_complete_coverage.txt)
num_genes_panel=$(wc -l < <(grep -v "chrY" ${MYOP_PANEL_BED}))
passed_coverage=false
# if [ "$num_genes_complete" -ge "$num_genes_panel" ]; then
if [ "$num_genes_to_complete" -eq 0 ]; then
    echo "✅ All genes in myopathy panel have complete coverage. Stop checking coverage."
    passed_coverage=true
else
    echo "Run coverage check for myopathy panel genes..."
    coverage_function
    num_genes_to_complete=$(wc -l < <(grep -v "chrY" ${TMP_PANEL}))
    echo "Genes with not enough coverage = ${num_genes_to_complete}"
    echo "Total genes in myopathy panel = ${num_genes_panel}"
    # if [ "$num_genes_complete" -ge "$num_genes_panel" ]; then
    if [ "$num_genes_to_complete" -eq 0 ]; then
        echo "✅ All genes in myopathy panel have complete coverage after this round."
        passed_coverage=true
    else
        echo "❌ Not all genes in myopathy panel have complete coverage. Continuing sequencing."
    fi
fi

# ---------------------------------------------------------------------------------------------
# LOD - Andre's script
# ---------------------------------------------------------------------------------------------
# Look through 4qA and see if the partia

count=$(awk -F'\t' '
NR==1 {
    # find column indexes
    for (i=1; i<=NF; i++) {
        if ($i == "MappedEstimatedCopies") c1 = i
        if ($i == "Haplotype") c2 = i
    }
    next
}

{
    gsub(/\r/, "", $0)  # remove Windows CR

    # count matches
    if ($c2 == "4qA" && $c1 > 3) {
        count++
    }
}

END {
    print count + 0   # ensure numeric output
}
' "$SUMMARY_TSV")

if [ "$count" -gt 0 ]; then
    echo "Found $count rows with 4qA and MappedEstimatedCopies > 3"

    lod=$(run_logged "spanning_ruleout" \
        "bash '${REALPATH}/spanning_ruleout_sample.sh' '$SUMMARY_TSV'")

    echo "LOD: $lod"
    lod_prob=$(echo $lod | cut -f2 -d" ")
    #
    lod_status=$(echo $lod | cut -f3 -d" ")
    if [ "$lod_status" == "TRUE" ]; then
        echo "✅ Enough sequencing to be $lod_prob confident of negative FSHD (≥ Threshold ($LOD_THRESH))"
    else
        echo "❌ LOD did not meet the threshold of $LOD_THRESH (Probability = $lod_prob)"
    fi
else
    echo "No 4qA rows with MappedEstimatedCopies > 3"
fi

# If all meet threshold, stop sequencing
if [[ "$passed_coverage" == "true" && \
      ( "$lod_status" == "TRUE" || "$negative_status" == "True" ) ]]; then
    echo "✅ Sample meets both coverage and negative/lod status criteria. STOP SEQUENCING."
    # exit 0
else
    echo "❌ Sample does not meet thresholds. Continuing sequencing."
    # exit 1
fi