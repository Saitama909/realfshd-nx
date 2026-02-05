#!/bin/bash

# INPUT_PATH=/data/neynjo24/ubam_paths/AS2603_bam_paths.tsv
# INPUT_PATH=/data/neynjo24/ubam_paths/ZU1108_bam_paths.tsv
INPUT_PATH=/data/neynjo24/ubam_paths/GUAT0705_bam_paths.tsv

INPUT_BLOW5=/data/neynjo24/blow5/QGXXOX250638_reads.blow5

source /data/neynjo24/realfshd_venv/bin/activate

counter=0
while IFS=$'\t' read -r SAMPLE_NAME OUTDIR INPUT_UBAM; do
    echo "Processing sample: $SAMPLE_NAME"
    echo "Output directory: $OUTDIR"
    echo "Input uBAM: $INPUT_UBAM"

    logs=/data/neynjo24/logs_realfshd/${SAMPLE_NAME}
    mkdir -p "$logs"

    chunk=$(basename $OUTDIR)
    logs_pre_file="${logs}/${SAMPLE_NAME}_${chunk}_fshd_run_prerun.log"

    time /data/neynjo24/realfshd/fshd/pre_mapping_n_rebasecall.sh "$INPUT_UBAM" "$OUTDIR" "$SAMPLE_NAME" "$INPUT_BLOW5" > ${logs_pre_file} 2>&1

    INPUT_UBAM_SUP="${OUTDIR}/${SAMPLE_NAME}_rebasecalled_supmeth_unaligned.bam"
    echo "Input uBAM sup: ${INPUT_UBAM_SUP}"
    if [ ! -f "$INPUT_UBAM_SUP" ]; then
        echo "Rebasecalled uBAM not found: ${INPUT_UBAM_SUP}"
        exit 1
    fi

    logs_file="${logs}/${SAMPLE_NAME}_${chunk}_fshd_run_script.log"
    time /data/neynjo24/realfshd/fshd/fshd_pipeline_wrapper.sh "$INPUT_UBAM_SUP" "$OUTDIR" "$SAMPLE_NAME" > ${logs_file} 2>&1

    # ----------------------------------------------------------------
    # Realtime stop analysis
    # 1. Merge myop panel BAM
    # 2. merge feature summary
    # 3. run fshd1 status generation
    # 4. run realtime_stop.sh
    PAR_DIR=$(dirname $OUTDIR)
    MERGED_MYOP_BAM="${PAR_DIR}/merged/${SAMPLE_NAME}_myop_panel.bam"
    CHUNK_MYOP_BAM="${OUTDIR}/${SAMPLE_NAME}_myop_panel.bam"

    MERGED_FEATURES_SUMMARY_TSV="${PAR_DIR}/merged/${SAMPLE_NAME}_mapped_features_summary.tsv"
    CHUNK_FEATURES_SUMMARY_TSV="${OUTDIR}/${SAMPLE_NAME}_mapped_features_summary.tsv"
    
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
extract_fshd1_status(summary_df, '${PAR_DIR}/merged/${SAMPLE_NAME}_fshd1_status_counts.tsv', '/data/neynjo24/realfshd/fshd/inputs/fshd1_status_template.tsv')
"

    time /data/neynjo24/realfshd/fshd/realtime_stop.sh \
        "${PAR_DIR}/merged/${SAMPLE_NAME}_fshd1_status_counts.tsv" \
        "${PAR_DIR}/merged/${SAMPLE_NAME}_mapped_features_summary.tsv" \
        "${PAR_DIR}/merged/${SAMPLE_NAME}_myop_panel.bam" > "${logs}/${SAMPLE_NAME}_${chunk}_fshd_realtime_stop.log" 2>&1
    
    # Copy all merged to merged_{counter} directory for record keeping
    OUTDIR_COUNTER="${PAR_DIR}/merged_${counter}"
    mkdir -p "${OUTDIR_COUNTER}"
    cp ${PAR_DIR}/merged/* ${OUTDIR_COUNTER}/.

    counter=$((counter + 1))

done < "$INPUT_PATH" #  <(head -n 6 "$INPUT_PATH")

deactivate