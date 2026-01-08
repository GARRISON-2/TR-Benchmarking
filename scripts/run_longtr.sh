#!/bin/bash
#SBATCH --job-name=longTR_one
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

module load nextflow
module load singularity/3.7.4

set -euo pipefail

# --- Args ---
# Usage: sbatch run_one_atarva.sbatch /abs/path/to/sample(.haplotagged).cram
CRAM="${1:?Usage: sbatch run_one_atarva.sbatch /path/to/sample.cram karyotype}"
KARYOTYPE="${2:?Usage: sbatch run_one_atarva.sbatch /path/to/sample.cram karyotype}"

# --- Config (match the Nextflow process inputs/paths) ---
REF="/pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta"
FAI="${REF}.fai"   # LongTR needs the FASTA; FAI is required by many tools; we validate below.
TR_REGIONS="/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/test-isolated-vc-catalog.longtr.bed"
LONGTR_BIN="/pl/active/dashnowlab/software/LongTR-1.2/LongTR"
OUTDIR="/pl/active/dashnowlab/projects/TR-benchmarking/longtr"

# Derived
mkdir -p "$OUTDIR" logs
base="$(basename "$CRAM")"
SAMPLE="${base%.haplotagged.cram}"
[[ "$SAMPLE" == "$base" ]] && SAMPLE="${base%.cram}"
OUT_VCF="${OUTDIR}/${SAMPLE}.longTR.vcf.gz"

echo "[$(date)] LongTR -> sample=$SAMPLE  cpus=${SLURM_CPUS_PER_TASK:-1}  out=$OUT_VCF"

# --- Environment ---
module load htslib/1.16
module load samtools/1.16 || true  # if available; useful to sanity-check CRAM/CRAI

MAX_TR_LEN="$(awk '{print $3-$2}' "$TR_REGIONS" | sort -n | tail -n 1)"
if [[ -z "${MAX_TR_LEN}" ]]; then
  echo "ERROR: Failed to compute MAX_TR_LEN from $TR_REGIONS" >&2
  exit 2
fi

# --- Alignment params from Nextflow (exact order preserved) ---
ALIGNMENT_PARAMS="-1.0,-0.458675,-1.0,-0.458675,-0.00005800168,-1,-1"

# --- Run LongTR (publishing to OUTDIR, like publishDir mode:'copy') ---
set -x
"$LONGTR_BIN" \
  --alignment-params "${ALIGNMENT_PARAMS}" \
  --fasta "${REF}" \
  --min-reads 1 \
  --max-tr-len "${MAX_TR_LEN}" \
  --regions "${TR_REGIONS}" \
  --bams "${CRAM}" \
  --tr-vcf "${OUT_VCF}"
set +x

echo "[$(date)] Done. Output: ${OUT_VCF}"
