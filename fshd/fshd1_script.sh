#!/bin/bash
#PBS -P kr68
#PBS -q normalbw
#PBS -l walltime=01:15:00
#PBS -l ncpus=3
#PBS -l mem=40GB
#PBS -l jobfs=10GB
#PBS -l wd

set -e
set -o pipefail

# Default inputs
INPUT_UBAM=""
INPUT_FASTQ=""

REALPATH=$(dirname "$(readlink -f "$0")")
REF=$REALPATH/genome/hs1.fa
echo "$REALPATH"
REGION_BED=${REALPATH}/inputs/d4z4_region.chm13.bed
FEATURES_FASTA=${REALPATH}/inputs/features.fasta
SHORT_FEATURES=${REALPATH}/inputs/short_features.fasta
PLAM=${REALPATH}/inputs/pLAM.fasta
PROBES=${REALPATH}/inputs/probes.fasta
PREFIX="SAMPLE"
OUTDIR="$PREFIX"
HAPLOTYPE_REFS=${REALPATH}/inputs/d4z4_repeats.fasta
MMI="${REALPATH}/inputs/hs1.mmi"
REMOVE_INTERMEDIATE_FILES=false
FSHD1_STATUS="${REALPATH}/inputs/fshd1_status_template.tsv"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --prefix) PREFIX="$2"; shift 2;;
        --outdir) OUTDIR="$2"; shift 2;;
        --input-ubam) INPUT_UBAM="$2"; shift 2;;
        --input-fastq) INPUT_FASTQ="$2"; shift 2;;
        --ref) REF="$2"; shift 2;;
        --region-bed) REGION_BED="$2"; shift 2;;
        --features-fasta) FEATURES_FASTA="$2"; shift 2;;
        --short-features-fasta) SHORT_FEATURES="$2"; shift 2;;
        --plam) PLAM="$2"; shift 2;;
        --haplotype-refs) HAPLOTYPE_REFS="$2"; shift 2;;
        --remove-intermediate-files) REMOVE_INTERMEDIATE_FILES="$2"; shift 2;;
        --minimap-ref-index) MMI="$2"; shift 2;;
        *) echo "Unknown option: $1"; exit 1;;
    esac
done

# Check if required arguments are set
if [[ (-z "$INPUT_UBAM" || ! -f "$INPUT_UBAM") && (-z "$INPUT_FASTQ" || ! -f "$INPUT_FASTQ") ]]; then
    echo "Error: Either --input-ubam or --input-fastq must be provided and must be a valid file."
    exit 1
fi

# Paths to conversion tools
PSLTOBED=pslToBed
BG2BW=bedGraphToBigWig

# Paths to tools
SAMTOOLS=samtools
MINIMAP2=minimap2
BEDTOOLS=bedtools

# Parameters
MAPQ=30
MIN_READ_LENGTH=30000
MODE=map-ont
# THREADS="${PBS_NCPUS:-8}"
THREADS=8
# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# -------------------- Added logging setup --------------------
LOG_DIR="/home/william/realfshd/pipeline_logs"
mkdir -p "$LOG_DIR"
PERF_LOG="${LOG_DIR}/fshd1_step_perf.log"

TEGRASTATS_BIN=$(command -v tegrastats || true)
GAWK_BIN=$(command -v gawk || true)

if [ ! -f "$PERF_LOG" ]; then
    echo -e "step\tstart_time\tend_time\twall_sec\tcpu_percent_time\ttegrastats_cpu_avg_percent\ttegrastats_cpu_max_percent\ttegrastats_gpu_avg_percent\ttegrastats_gpu_max_percent\tstatus" | tee "$PERF_LOG"
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
        "$TEGRASTATS_BIN" --interval 200 > "$TEGRA_TMP" 2>/dev/null &
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

echo "Convert ubam/fastq to fastq/ubam"
ubam_exist=false
fastq_exist=false
if [[ (-n "$INPUT_UBAM") && (-f "$INPUT_UBAM") ]]; then
    ubam_exist=true
    UBAM_BASE=$(basename "$INPUT_UBAM")
    run_logged "samtools_fastq" \
        "$SAMTOOLS fastq -TMM,ML \"$INPUT_UBAM\" > \"${OUTDIR}/${UBAM_BASE%.bam}.fastq\""
    INPUT_FASTQ="${OUTDIR}/${UBAM_BASE%.bam}.fastq"
elif [[ (-n "$INPUT_FASTQ") && (-f "$INPUT_FASTQ") ]]; then
    fastq_exist=true
    FASTQ_BASE=$(basename "$INPUT_FASTQ")
    run_logged "samtools_import" \
        "$SAMTOOLS import -u -s \"$INPUT_FASTQ\" -o \"${OUTDIR}/${FASTQ_BASE%.fastq}.bam\""
    INPUT_UBAM="${OUTDIR}/${FASTQ_BASE%.fastq}.bam"
fi

echo "MERYL"
REF_NAME=$(basename "$REF" .fa)
MERYL_DB="${REALPATH}/inputs/merylDB_${REF_NAME}"
REPETITIVE_REGIONS="${REALPATH}/inputs/${REF_NAME}_repetitive_k15.txt"

if [ ! -d "$MERYL_DB" ]; then
    echo "Creating meryl database..."
    run_logged "meryl_count" \
        "meryl count k=15 output \"$MERYL_DB\" \"$REF\""
else
    echo "Meryl database already exists at $MERYL_DB. Skipping this step."
fi

if [ ! -f "$REPETITIVE_REGIONS" ]; then
    echo "Extracting repetitive regions..."
    run_logged "meryl_print" \
        "meryl print greater-than distinct=0.9998 \"$MERYL_DB\" > \"$REPETITIVE_REGIONS\""
else
    echo "Repetitive regions file already exists at $REPETITIVE_REGIONS. Skipping this step."
fi

echo "$(realpath "$OUTDIR")"
echo "$PREFIX"

echo "Extracting myopathy panel from fast basecall chm13 mapping"
MYOP_PANEL=${REALPATH}/inputs/myopathy_panel_draft_2.gene_targets.chm13.bed
run_logged "extract_myop_panel" \
    "$SAMTOOLS view -L \"$MYOP_PANEL\" -Sb \"${OUTDIR}/../${PREFIX}_fast_chm13.bam\" | $SAMTOOLS sort -o \"${OUTDIR}/${PREFIX}_myop_panel.bam\""
run_logged "index_myop_panel" \
    "$SAMTOOLS index \"${OUTDIR}/${PREFIX}_myop_panel.bam\""

echo "minimap winnowmap"
SPLIT_PREFIX="${OUTDIR}/${PREFIX}_hs1_split_tmp"

run_logged "minimap2_chm13" \
    "$MINIMAP2 -x \"$MODE\" -a \"$MMI\" \"$INPUT_FASTQ\" --split-prefix \"$SPLIT_PREFIX\" -t ${THREADS} | $SAMTOOLS sort -o \"${OUTDIR}/${PREFIX}_chm13.bam\""

run_logged "extract_d4z4_reads" \
    "$SAMTOOLS view -L \"$REGION_BED\" -Sb \"${OUTDIR}/${PREFIX}_chm13.bam\" | $SAMTOOLS sort -o \"${OUTDIR}/${PREFIX}_d4z4_reads_of_interest.bam\""
run_logged "index_d4z4_reads" \
    "$SAMTOOLS index \"${OUTDIR}/${PREFIX}_d4z4_reads_of_interest.bam\""

READS_OF_INTEREST_BAM="${OUTDIR}/${PREFIX}_reads_of_interest.bam"

# if [[ -s "$READS_OF_INTEREST_BAM" && -s "${READS_OF_INTEREST_BAM}.bai" ]]; then
#     echo "Using existing reads_of_interest BAM: ${READS_OF_INTEREST_BAM}"
# else
#     echo "Remap with winnowmap to CHM13"
#     run_logged "winnowmap_reads_of_interest" \
#         "$SAMTOOLS view -N <($SAMTOOLS view \"${OUTDIR}/${PREFIX}_d4z4_reads_of_interest.bam\" | cut -f1 | sort | uniq) -b \"${INPUT_UBAM}\" | \
#          $SAMTOOLS fastq -TMM,ML - | \
#          winnowmap -W \"${REPETITIVE_REGIONS}\" -Y -y -ax \"$MODE\" \"$REF\" - -t ${THREADS} | \
#          $SAMTOOLS view -L \"$REGION_BED\" -Sb | \
#          $SAMTOOLS sort -o \"${READS_OF_INTEREST_BAM}\""
#     run_logged "index_reads_of_interest" \
#         "$SAMTOOLS index \"${READS_OF_INTEREST_BAM}\""
# fi
echo "Remap with winnowmap to CHM13"
run_logged "winnowmap_reads_of_interest" \
    "$SAMTOOLS view -N <($SAMTOOLS view \"${OUTDIR}/${PREFIX}_d4z4_reads_of_interest.bam\" | cut -f1 | sort | uniq) -b \"${INPUT_UBAM}\" | \
        $SAMTOOLS fastq -TMM,ML - | \
        winnowmap -W \"${REPETITIVE_REGIONS}\" -Y -y -ax \"$MODE\" \"$REF\" - -t ${THREADS} | \
        $SAMTOOLS view -L \"$REGION_BED\" -Sb | \
        $SAMTOOLS sort -o \"${READS_OF_INTEREST_BAM}\""
run_logged "index_reads_of_interest" \
    "$SAMTOOLS index \"${READS_OF_INTEREST_BAM}\""
    
echo "Generating FASTA from uBAM"
run_logged "generate_reads_fasta" \
    "$SAMTOOLS view -N <($SAMTOOLS view \"${OUTDIR}/${PREFIX}_reads_of_interest.bam\" | cut -f1 | sort | uniq) -b \"$INPUT_UBAM\" | $SAMTOOLS fasta > \"${OUTDIR}/${PREFIX}_reads_of_interest.fasta\""

echo "Converting BAM to BED"
run_logged "bam_to_bed" \
    "bedtools bamtobed -i \"${OUTDIR}/${PREFIX}_reads_of_interest.bam\" > \"${OUTDIR}/${PREFIX}_reads_of_interest.bed\""

run_logged "minimap_features" \
    "${REALPATH}/minimap_features.sh \"$OUTDIR\" \"$PREFIX\" \"${OUTDIR}/${PREFIX}_reads_of_interest.fasta\" \"$FEATURES_FASTA\" \"$SHORT_FEATURES\" \"$PLAM\""

run_logged "seqkit_sslp" \
    "seqkit amplicon --bed -F GGTGGAGTTCTGGTTTCAGC -R CCTGTGCTTCAGAGGCATTTG -m 2 \"${OUTDIR}/${PREFIX}_reads_of_interest.fasta\" > \"${OUTDIR}/${PREFIX}_SSLP.bed\""

echo "Aligning reads classified to specific haplotypes"
aligned_bams=()

for chrom in chr4 chr10; do
    if [[ "$chrom" == "chr4" ]]; then
        haplotype="4A"
    else
        haplotype="10A"
    fi

    haplotype_ref="${OUTDIR}/${haplotype}_reference.fasta"
    haplotype_reads="${OUTDIR}/${PREFIX}_${haplotype}_reads.fastq"
    aligned_bam="${OUTDIR}/${PREFIX}_${haplotype}_aligned.bam"

    run_logged "extract_${haplotype}_reference" \
        "grep -A1 \"${haplotype}\" \"$HAPLOTYPE_REFS\" > \"$haplotype_ref\""

    run_logged "extract_${haplotype}_reads" \
        "seqtk subseq \"$INPUT_FASTQ\" <(grep ${chrom} \"${OUTDIR}/${PREFIX}_reads_of_interest.bed\" | cut -f4 | sort | uniq) > \"$haplotype_reads\""

    run_logged "align_${haplotype}_reads" \
        "$MINIMAP2 --secondary=no --MD -ax map-ont \"$haplotype_ref\" \"$haplotype_reads\" | $SAMTOOLS sort -o \"$aligned_bam\""

    run_logged "index_${haplotype}_bam" \
        "$SAMTOOLS index \"$aligned_bam\""

    aligned_bams+=("$aligned_bam")
    rm -f "$haplotype_ref" "$haplotype_reads"
done
echo "Haplotype-specific alignment complete."

echo "Merging haplotype-aligned BAM files"
merged_bam="${OUTDIR}/${PREFIX}_aligned_haplotypes.bam"
run_logged "merge_haplotype_bams" \
    "$SAMTOOLS merge -f \"$merged_bam\" ${aligned_bams[*]}"
run_logged "index_merged_haplotype_bam" \
    "$SAMTOOLS index \"$merged_bam\""

rm "${aligned_bams[@]}" "${aligned_bams[@]/%.bam/.bam.bai}"

echo "Parsing PSL to generate summary table"
run_logged "read_classification" \
    "python3 \"$REALPATH/helper/read_classification.py\" \
    --psl \"${OUTDIR}/${PREFIX}_mapped_features.psl\" \
    --bed \"${OUTDIR}/${PREFIX}_reads_of_interest.bed\" \
    --output \"${OUTDIR}/${PREFIX}_mapped_features_summary.tsv\" \
    --fasta \"${OUTDIR}/${PREFIX}_reads_of_interest.fasta\" \
    --sslp \"${OUTDIR}/${PREFIX}_SSLP.bed\" \
    --aligned-bam \"${OUTDIR}/${PREFIX}_aligned_haplotypes.bam\""

echo "Generating ordered alignment sequences"
alignment_script="$REALPATH/helper/alignment_order_clipping.py"
run_logged "alignment_order_clipping" \
    "python3 \"$alignment_script\" \
    --bam \"${OUTDIR}/${PREFIX}_aligned_haplotypes.bam\" \
    --fasta \"${OUTDIR}/${PREFIX}_reads_of_interest.fasta\" \
    --output_table \"${OUTDIR}/${PREFIX}_d4z4_units.tsv\" \
    --xapi_bed \"${OUTDIR}/${PREFIX}_xapi_sites.bed\" \
    --blni_bed \"${OUTDIR}/${PREFIX}_blni_sites.bed\" \
    --main_tsv \"${OUTDIR}/${PREFIX}_mapped_features_summary.tsv\" \
    --repeats_bed \"${OUTDIR}/${PREFIX}_d4z4_repeats.bed\""

run_logged "potential_hybrid_awk" \
    "$SAMTOOLS view \"${OUTDIR}/${PREFIX}_reads_of_interest.bam\" | awk '{
        primary_chr = \$3
        for (i=12; i<=NF; i++) {
            if (\$i ~ /^SA:Z:/) {
                sub(/^SA:Z:/, \"\", \$i)
                n = split(\$i, sa_alignments, \";\")
                for (j=1; j<=n; j++) {
                    if (sa_alignments[j] == \"\") continue
                    split(sa_alignments[j], sa_fields, \",\")
                    sa_chr = sa_fields[1]
                    mapq = sa_fields[5]
                    if (primary_chr != sa_chr && mapq >= 30) {
                        print \$1, primary_chr, sa_alignments[j]
                    }
                }
            }
        }
    }' > \"${OUTDIR}/${PREFIX}_potential_hybrid.tsv\""

run_logged "flag_repeats" \
    "python3 \"$REALPATH/helper/flag_repeats.py\" \
    --main_tsv \"${OUTDIR}/${PREFIX}_mapped_features_summary.tsv\" \
    --repeats_bed \"${OUTDIR}/${PREFIX}_d4z4_repeats.bed\" \
    --hybrid \"${OUTDIR}/${PREFIX}_potential_hybrid.tsv\" \
    --fasta \"${OUTDIR}/${PREFIX}_reads_of_interest.fasta\""

run_logged "add_colors_to_bed" \
    "python3 \"$REALPATH/helper/add_colors_to_bed.py\" \
    --main_tsv \"${OUTDIR}/${PREFIX}_mapped_features_summary.tsv\" \
    --features_bed \"${OUTDIR}/${PREFIX}_mapped_features.bed\" \
    --repeats_bed \"${OUTDIR}/${PREFIX}_d4z4_repeats.bed\" \
    --sslp_bed \"${OUTDIR}/${PREFIX}_SSLP.bed\" \
    --output_bed \"${OUTDIR}/${PREFIX}_all_features.bed\""

run_logged "distal_haplotype_blast" \
    "${REALPATH}/distal_haplotype_blast.sh \"${OUTDIR}\" \"${PREFIX}\""

run_logged "generate_meth_reads_bam" \
    "$SAMTOOLS view -N <($SAMTOOLS view \"${OUTDIR}/${PREFIX}_reads_of_interest.bam\" | cut -f1 | sort | uniq) -b \"${INPUT_UBAM}\" | \
     $SAMTOOLS fastq -TMM,ML - | \
     $MINIMAP2 -t ${THREADS} -Y -y -x asm5 -a --secondary=no \"${OUTDIR}/${PREFIX}_reads_of_interest.fasta\" - | \
     $SAMTOOLS view -F 2308 -b - | \
     $SAMTOOLS sort -@ ${THREADS} - > \"${OUTDIR}/${PREFIX}_meth_reads.bam\""

run_logged "index_meth_reads_bam" \
    "$SAMTOOLS index \"${OUTDIR}/${PREFIX}_meth_reads.bam\""

run_logged "minimod_view" \
    "minimod view -c m \"${OUTDIR}/${PREFIX}_reads_of_interest.fasta\" \"${OUTDIR}/${PREFIX}_meth_reads.bam\" | awk 'NR > 1 {print \$4\"\t\"\$5\"\t\"\$5+1\"\t.\t\"\$7}' > \"${OUTDIR}/${PREFIX}_meth_reads.bed\""

fshd1_read_status_tsv="${OUTDIR}/${PREFIX}_fshd1_status_counts.tsv"
run_logged "methylation_summary" \
    "python3 \"$REALPATH/helper/methylation_summary.py\" \
    --main_tsv \"${OUTDIR}/${PREFIX}_mapped_features_summary.tsv\" \
    --features_bed \"${OUTDIR}/${PREFIX}_all_features.bed\" \
    --meth \"${OUTDIR}/${PREFIX}_meth_reads.bed\" \
    --output \"${OUTDIR}/${PREFIX}_methylation_summary.tsv\" \
    --fshd1_read_status_tsv \"$fshd1_read_status_tsv\" \
    --status_template_tsv \"$FSHD1_STATUS\" \
    --updated_bed \"${OUTDIR}/${PREFIX}_updated_features.bed\""

run_logged "replace_summary_tsv" \
    "mv \"${OUTDIR}/${PREFIX}_methylation_summary.tsv\" \"${OUTDIR}/${PREFIX}_mapped_features_summary.tsv\""

if [ ! -f "${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes" ]; then
    echo "Creating ${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes in $OUTDIR..."
    if [ ! -f "${OUTDIR}/${PREFIX}_reads_of_interest.fasta.fai" ]; then
        run_logged "faidx_reads_fasta" \
            "$SAMTOOLS faidx \"${OUTDIR}/${PREFIX}_reads_of_interest.fasta\""
    fi
    run_logged "create_chrom_sizes" \
        "cut -f1,2 \"${OUTDIR}/${PREFIX}_reads_of_interest.fasta.fai\" > \"${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes\""
else
    echo "${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes already exists in $OUTDIR."
fi

run_logged "high_prob_bedgraph" \
    "awk '\$5 > 0.75 {print \$1, \$2, \$3, \$5}' \"${OUTDIR}/${PREFIX}_meth_reads.bed\" | sort -k1,1 -k2,2n > \"${OUTDIR}/${PREFIX}_meth_reads.high_prob.bedGraph\""

run_logged "high_prob_bigwig" \
    "$BG2BW \"${OUTDIR}/${PREFIX}_meth_reads.high_prob.bedGraph\" \"${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes\" \"${OUTDIR}/${PREFIX}_meth_reads.high_prob.bw\""

run_logged "low_prob_bedgraph" \
    "awk '\$5 < 0.25 {print \$1, \$2, \$3, \$5}' \"${OUTDIR}/${PREFIX}_meth_reads.bed\" | sort -k1,1 -k2,2n > \"${OUTDIR}/${PREFIX}_meth_reads.low_prob.bedGraph\""

run_logged "low_prob_bigwig" \
    "$BG2BW \"${OUTDIR}/${PREFIX}_meth_reads.low_prob.bedGraph\" \"${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes\" \"${OUTDIR}/${PREFIX}_meth_reads.low_prob.bw\""

run_logged "undefined_prob_bedgraph" \
    "awk '\$5 >= 0.25 && \$5 <= 0.75 {print \$1, \$2, \$3, \$5}' \"${OUTDIR}/${PREFIX}_meth_reads.bed\" | sort -k1,1 -k2,2n > \"${OUTDIR}/${PREFIX}_meth_reads.undefined.bedGraph\""

run_logged "undefined_prob_bigwig" \
    "$BG2BW \"${OUTDIR}/${PREFIX}_meth_reads.undefined.bedGraph\" \"${OUTDIR}/${PREFIX}_reads_of_interest.chrom.sizes\" \"${OUTDIR}/${PREFIX}_meth_reads.undefined.bw\""

run_logged "cleanup_bedgraphs" \
    "rm \"${OUTDIR}/${PREFIX}_meth_reads.high_prob.bedGraph\" \"${OUTDIR}/${PREFIX}_meth_reads.low_prob.bedGraph\" \"${OUTDIR}/${PREFIX}_meth_reads.undefined.bedGraph\""

if [[ "$REMOVE_INTERMEDIATE_FILES" == "true" ]]; then
    run_logged "cleanup_intermediate" \
        "rm ${OUTDIR}/${PREFIX}_*blast* ${OUTDIR}/${PREFIX}_*top* ${OUTDIR}/${PREFIX}_*psl ${OUTDIR}/${PREFIX}_*d4z4*"
fi

if [[ "$ubam_exist" == "true" ]]; then
    run_logged "cleanup_temp_fastq" \
        "rm \"${INPUT_FASTQ}\""
elif [[ "$fastq_exist" == "true" ]]; then
    run_logged "cleanup_temp_ubam" \
        "rm \"${INPUT_UBAM}\""
fi