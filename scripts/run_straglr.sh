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

# --- Args ---
# Usage: sbatch run_one_atarva.sbatch /abs/path/to/sample(.haplotagged).cram
CRAM="${1:?Usage: sbatch run_one_atarva.sbatch /path/to/sample.cram karyotype}"
KARYOTYPE="${2:?Usage: sbatch run_one_atarva.sbatch /path/to/sample.cram karyotype}"

# --- Config ---
REF="/pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta"
BED="/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/STRchive-disease-loci.hg38.longTR.bed"
OUTDIR="/pl/active/dashnowlab/projects/TR-benchmarking/straglr"


mkdir -p "$OUTDIR" logs

# --- Environment ---
module load miniforge/24.11.3-0
conda activate straglr_1.5.5

# --- Derive sample name ---
base="$(basename "$CRAM")"
SAMPLE="${base%.haplotagged.cram}"
[[ "$SAMPLE" == "$base" ]] && SAMPLE="${base%.cram}"

OUT_VCF="${OUTDIR}/${SAMPLE}.atarva-3.0.1.vcf"

echo "[$(date)] straglr -> sample=$SAMPLE  threads=4"
 # 2) Run Straglr (targeted loci genotyping)
python3 /pl/active/dashnowlab/software/straglr/straglr.py \
  "$CRAM" \
  "$REF" \
  "${OUTDIR}/${SAMPLE}" \
  --loci "${BED}" \
  --chroms 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y \
  --min_support 1 \
  --min_cluster_size 1 \
  --max_num_clusters 2 \
  --genotype_in_size \
  --nprocs 4