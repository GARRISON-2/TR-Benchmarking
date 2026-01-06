#!/bin/bash
#SBATCH --job-name=atarva_one
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

#This is how to run it
#while IFS= read -r CRAM; do   [[ -z "$CRAM" || "$CRAM" =~ ^[[:space:]]*# ]] && continue;   SAMPLE="$(basename "$CRAM" .haplotagged.cram)";   SAMPLE="${SAMPLE%.cram}";   sbatch -J "atarva_${SAMPLE}" 2_run_atarva_sbatch.sh "$CRAM" "$KARYOTYPE"; done < haplotagged_test.list

# --- Args ---
# Usage: sbatch run_one_atarva.sbatch /abs/path/to/sample(.haplotagged).cram
CRAM="${1:?Usage: sbatch run_one_atarva.sbatch /path/to/sample.cram karyotype}"
KARYOTYPE="${2:?Usage: sbatch run_one_atarva.sbatch /path/to/sample.cram karyotype}"

# --- Config ---
REF="/pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta"
#BED="catalogs/STRchive-disease-loci.hg38.atarva.bed.gz"
BED="/pl/active/dashnowlab/work/ealiyev/software/trxplorer_catalogs/TR_catalog.40885_62_loci.atarva.bed.gz"
OUTDIR="/pl/active/dashnowlab/projects/TR-benchmarking/atarva"
THREADS="${SLURM_CPUS_PER_TASK:-8}"


mkdir -p "$OUTDIR" logs

# --- Environment ---
module load bcftools/1.16
module load bedtools/2.29.1
module load htslib/1.16
module load miniforge/24.11.3-0
conda activate atarva-3.0.1

# --- Derive sample name ---
base="$(basename "$CRAM")"
SAMPLE="${base%.haplotagged.cram}"
[[ "$SAMPLE" == "$base" ]] && SAMPLE="${base%.cram}"

OUT_VCF="${OUTDIR}/${SAMPLE}.atarva-3.0.1.vcf"

echo "[$(date)] atarva -> sample=$SAMPLE  threads=$THREADS"
atarva -t "$THREADS" \
 -f "$REF" \
 -b "$CRAM" \
 -r "$BED" \
 --format cram \
 --haplotag HP \
 --min-reads 1 \
 --decompose \
 --karyotype $KARYOTYPE \
 -o "$OUT_VCF"