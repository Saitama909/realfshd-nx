#!/bin/bash

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

# Convert to FASTQ
UBAM_BASE=$(basename $INPUT_UBAM)
INPUT_FASTQ="${OUTDIR}/${UBAM_BASE%.bam}.fastq"
$SAMTOOLS fastq -TMM,ML "$INPUT_UBAM" > "${INPUT_FASTQ}"

echo "minimap"
CHM13_BAM="${OUTDIR}/${PREFIX}_fast_chm13.bam"
$MINIMAP2 -x $MODE -a "$MMI" "$INPUT_FASTQ" -t ${PBS_NCPUS:-8} |\
    $SAMTOOLS sort -o "${CHM13_BAM}"
$SAMTOOLS index "${CHM13_BAM}"

echo "extracting reads mapping to the region of interest"
# Get d4z4 reads
D4Z4_ROI_BAM="${OUTDIR}/${PREFIX}_d4z4_reads_of_interest_fast.bam"
$SAMTOOLS view -L "$REGION_BED" -Sb "${CHM13_BAM}" | \
    $SAMTOOLS sort -o "${D4Z4_ROI_BAM}"
$SAMTOOLS index "${D4Z4_ROI_BAM}"

# -----------------------------------------------------------------------------------
#  Extract D4Z4 reads for rebasecalling with methylation
# -----------------------------------------------------------------------------------
d4z4_rid="${OUTDIR}/${PREFIX}_d4z4_rid.txt"

$SAMTOOLS view "${D4Z4_ROI_BAM}" | cut -f1 | sort -u > "$d4z4_rid"

# Extract reads from uBAM
D4Z4_UBAM="${OUTDIR}/${PREFIX}_d4z4_reads_in_ubam.bam"
echo "extract d4z4 reads from ubam"
$SAMTOOLS view -N "$d4z4_rid" -b "$INPUT_UBAM" > "${D4Z4_UBAM}"
$SAMTOOLS index "${D4Z4_UBAM}"

# Extract read IDs with and without pi tags from uBAM
# d4z4_rid="${OUTDIR}/${PREFIX}_d4z4_rid.txt"

# Extract reads without pi tags
grep -vP 'pi:Z:\K\S+' <( $SAMTOOLS view "${D4Z4_UBAM}" ) | cut -f1 | sort -u > "$d4z4_rid"

# # Extract read with pi tags
grep -oP 'pi:Z:\K\S+' <( $SAMTOOLS view "${D4Z4_UBAM}" ) | sort -u >> "$d4z4_rid"

sort -u "$d4z4_rid" > "${OUTDIR}/d4z4_rid_final.txt"
mv "${OUTDIR}/d4z4_rid_final.txt" "$d4z4_rid"

# Index INPUT BLOW5
${SLOW5TOOLS} index ${INPUT_BLOW5}

# After filtering reads with on D4Z4
echo "extract d4z4 reads from blow5"
EXTRACTED_BLOW5="${OUTDIR}/${PREFIX}_d4z4_reads.blow5"
${SLOW5TOOLS} get "${INPUT_BLOW5}" --list "${d4z4_rid}" -o "${EXTRACTED_BLOW5}"
${SLOW5TOOLS} index ${EXTRACTED_BLOW5}

# # rebasecall with SUP
echo "rebasecall with SUP"
${SLOW5_DORADO} basecaller -x cuda:all --models-directory ${DORADO_MODEL_PATH} ${DORADO_MODEL_PATH}/${DORADO_MODEL_SUP} "${EXTRACTED_BLOW5}" \
    --modified-bases 5mCG_5hmCG \
    > ${REBASECALLED_UBAM}

# Remove blow5
rm "${EXTRACTED_BLOW5}" "${EXTRACTED_BLOW5}.idx"
rm "${INPUT_FASTQ}"
rm "${D4Z4_ROI_BAM}" "${D4Z4_ROI_BAM}.bai"
# rm "$d4z4_rid"