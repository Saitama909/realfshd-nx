#!/bin/bash

# ----------------------------------------------------------------
# Realtime stop analysis
# 1. Merge myop panel BAM
# 2. merge feature summary
# 3. run fshd1 status generation
# 4. run realtime_stop.sh
OUTDIR=$1
PAR_DIR=$(dirname $OUTDIR)
SAMPLE_NAME=SAMPLE
MERGED_MYOP_BAM="${PAR_DIR}/merged/${SAMPLE_NAME}_myop_panel.bam"
CHUNK_MYOP_BAM="${OUTDIR}/${SAMPLE_NAME}_myop_panel.bam"

MERGED_FEATURES_SUMMARY_TSV="${PAR_DIR}/merged/${SAMPLE_NAME}_mapped_features_summary.tsv"
CHUNK_FEATURES_SUMMARY_TSV="${OUTDIR}/${SAMPLE_NAME}_mapped_features_summary.tsv"

REALPATH=$(dirname "$(readlink -f "$0")")

echo "Merging myopathy panel BAM to: ${MERGED_MYOP_BAM}"
echo "Merging features summary to: ${MERGED_FEATURES_SUMMARY_TSV}"

if [ ! -d "${PAR_DIR}/merged" ]; then
    mkdir -p "${PAR_DIR}/merged"
    cp "${CHUNK_MYOP_BAM}" "${MERGED_MYOP_BAM}"
    samtools index "${MERGED_MYOP_BAM}"
    cp "${CHUNK_FEATURES_SUMMARY_TSV}" "${MERGED_FEATURES_SUMMARY_TSV}"
else
    # Append to existing merged summary file (excluding header)
    tail -n +2 "${CHUNK_FEATURES_SUMMARY_TSV}" >> "${MERGED_FEATURES_SUMMARY_TSV}"

    # Merge BAM files
    tmp_bam="${MERGED_MYOP_BAM%.bam}.tmp.bam"
    samtools merge -f "$tmp_bam" \
        "$MERGED_MYOP_BAM" \
        "$CHUNK_MYOP_BAM"
    mv "$tmp_bam" "$MERGED_MYOP_BAM"
    # samtools merge -f "${MERGED_MYOP_BAM}" "${MERGED_MYOP_BAM}" "${CHUNK_MYOP_BAM}"
    samtools index "${MERGED_MYOP_BAM}"
fi

# FSHD1 status generation
python3 -c "
from helper.methylation_summary import extract_fshd1_status
import pandas as pd
summary_df = pd.read_csv('${MERGED_FEATURES_SUMMARY_TSV}', sep='\t')
extract_fshd1_status(summary_df, '${PAR_DIR}/merged/${SAMPLE_NAME}_fshd1_status_counts.tsv', '${REALPATH}/inputs/fshd1_status_template.tsv')
"

time ${REALPATH}/realtime_stop.sh \
    "${PAR_DIR}/merged/${SAMPLE_NAME}_fshd1_status_counts.tsv" \
    "${PAR_DIR}/merged/${SAMPLE_NAME}_mapped_features_summary.tsv" \
    "${PAR_DIR}/merged/${SAMPLE_NAME}_myop_panel.bam" > "${logs}/${SAMPLE_NAME}_${chunk}_fshd_realtime_stop.log" 2>&1

# Copy all merged to merged_{counter} directory for record keeping
OUTDIR_COUNTER="${PAR_DIR}/merged_${counter}"
mkdir -p "${OUTDIR_COUNTER}"
cp ${PAR_DIR}/merged/* ${OUTDIR_COUNTER}/.



