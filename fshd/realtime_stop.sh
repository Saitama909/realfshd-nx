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
# Set thresholds
# ---------------------------------------------------------------------------------------------
POSITIVE_THRESH=3 # Positive FSHD
NEGATIVE_THRESH=5 # Negative FSHD
COVERAGE_THRESH=30 # Threshold for gene to stop being in read fish from panel
LOD_THRESH=0.95 # Limit of detection threshold for FSHD

# ---------------------------------------------------------------------------------------------
# Positive case - Stop sequencing
# ---------------------------------------------------------------------------------------------
awk -v thresh="$POSITIVE_THRESH" -F'\t' '
$1=="Diagnostic/Positive" && $7>=thresh {
    print "✅ Positive FSHD1 case, enough sequencing"
    exit 100
}' "${STATUS_TSV}"

# Check awk exit code
# Stop the script if it's enough
if [ $? -eq 100 ]; then
    exit 0   # stop the script
fi

# ---------------------------------------------------------------------------------------------
# Easy negative case (1x 4qA with more than 10 copies and 1x 4qB)
# ---------------------------------------------------------------------------------------------
read negative_status message < <(
awk -F'\t' '
$1 == "Diagnostic/Negative" { 
    count++; 
    if ($NF >= 3) ok++ 
} 
END { 
    if (count == 2 && ok == 2) {
        print "True", "✅ Have a non-permissive 4qA and 4qB reads with ≥3 reads each"; 
    } else {
        print "False", "❌ One or both Diagnostic/Negative rows have <3 reads" 
    }
}' $STATUS_TSV
)

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
    samtools depth -a -b $TMP_PANEL $MYOP_BAM > $COVERAGE_OUTPUT

    # Currently col 2 is coordinate and col 3 is depth
    # want to make it col 2 and col 3 is coordinate (duplicated col2)
    BED_COV="${SAMPLE_DIR}/myop_panel_chm13_coverage.bed"
    awk -F'\t' '{print $1"\t"$2"\t"$2"\t"$3}' $COVERAGE_OUTPUT > "$BED_COV"

    echo "INTERSECT"
    # get intersected regions of coverage and panel
    bedtools intersect -a $TMP_PANEL -b $BED_COV -wa -wb > ${SAMPLE_DIR}/panel_cov_intersect.bed

    awk -v cov_thresh="$COVERAGE_THRESH" \
        -F'\t' '
    {
        gene = $4
        cov  = $8

        sum_cov[gene] += cov
        base_count[gene]++

        if (cov < cov_thresh)
            under[gene]++
    }

    END {
        printf "%-20s\t%-15s\t%-10s\t%-15s\t%-18s\t%-15s\t%s\n",
            "Gene", "Avg_Coverage", "Length", "Total_Coverage",
            "Bases_Under_Threshold", "Pct_Bases_Above", "Pass_Threshold"

        for (g in sum_cov) {
            total = base_count[g]
            under_count = (g in under ? under[g] : 0)
            above_pct = 1 - (under_count / total)
            avg = sum_cov[g] / total

            pass = (avg >= cov_thresh && above_pct >= 0.95 ? "YES" : "NO")

            printf "%-20s\t%-15.2f\t%-10d\t%-15d\t%-18d\t%-15.3f\t%s\n",
                g, avg, total, sum_cov[g], under_count, above_pct, pass
        }
    }' "${SAMPLE_DIR}/panel_cov_intersect.bed" > "$output"

    # Remove genes with complete coverage from tmp panel bed file
    while IFS=$'\t' read -r gene cov_avg length  total_cov    bases_under_thresh prop_above_thresh   pass_thresh;
    do
        # echo "Removing gene with complete coverage: $gene"
        # echo "Gene: $gene, Pass Threshold: $pass_thresh"
        if [ "$pass_thresh" == "YES" ]; then
            # gene="$(echo "$gene" | tr -d '[:space:]')"
            gene="$(echo "$gene" | sed -e 's/[[:space:]]*$//')"
            # echo "Gene $gene has complete coverage. Removing from temp panel."
            grep -v "$gene" $TMP_PANEL > "${TMP_PANEL}.tmp"
            mv "${TMP_PANEL}.tmp" $TMP_PANEL
        fi
    done < ${SAMPLE_DIR}/genes_complete_coverage.txt
}
    # if (print_header == 1) {
    #     }

# if lines in tmp panel = lines in MYOP PANEL
# then all genes have complete coverage
if [ -f ${TMP_PANEL} ]; then
    num_genes_to_complete=$(wc -l < ${TMP_PANEL})
else
    num_genes_to_complete=$(wc -l < $MYOP_PANEL_BED)
fi
# if [ -f ${SAMPLE_DIR}/genes_complete_coverage.txt ]; then
#     num_genes_complete=$(wc -l < ${SAMPLE_DIR}/genes_complete_coverage.txt)
# else
#     num_genes_complete=0
# fi
# num_genes_complete=$(wc -l < ${SAMPLE_DIR}/genes_complete_coverage.txt)
num_genes_panel=$(wc -l < $MYOP_PANEL_BED)
passed_coverage=false
# if [ "$num_genes_complete" -ge "$num_genes_panel" ]; then
if [ "$num_genes_to_complete" -eq 0 ]; then
    echo "✅ All genes in myopathy panel have complete coverage. Stop checking coverage."
    passed_coverage=true
else
    echo "Run coverage check for myopathy panel genes..."
    coverage_function
    num_genes_to_complete=$(wc -l < ${TMP_PANEL})
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
    lod=$(bash ${REALPATH}/spanning_ruleout_sample.sh "$SUMMARY_TSV")

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
    echo "✅ Sample meets both coverage and negative/lod status criteria. Stopping sequencing."
    # exit 0
else
    echo "❌ Sample does not meet thresholds. Continuing sequencing."
    # exit 1
fi
