#!/bin/bash
#SBATCH --job-name=trsv_call
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

# --- Input arguments ---
CRAM="${1:?Usage: sbatch run_trsv.sbatch /path/to/sample.cram}"
REF="/pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta"     # change if needed
PLATFORM="ont"                              # ont | hifi | clr
BUILD="38"                                  # 37 | 38 | T2T
THREADS=8

# --- Output prefix (derived from CRAM filename) ---
SAMPLE=$(basename "$CRAM" .haplotagged.cram)
OUT_PREFIX="/pl/active/dashnowlab/projects/TR-benchmarking/trsv/${SAMPLE}_trsv"

# --- Paths ---
SIF="/projects/ealiyev@xsede.org/TRsv-v1.1.sif"

# --- Load Singularity module ---
module load singularity || true

# --- Bind directories ---
DATA_DIR=$(dirname "$CRAM")
REF_DIR=$(dirname "$REF")
OUT_DIR=$(dirname "$OUT_PREFIX")

echo "=============================="
echo "[TRSV JOB INFO]"
echo "Input CRAM:     $CRAM"
echo "Reference:      $REF"
echo "Output prefix:  $OUT_PREFIX"
echo "Sample name:    $SAMPLE"
echo
echo "[BIND PATHS]"
echo "Data dir:       $DATA_DIR  -->  /data"
echo "Reference dir:  $REF_DIR   -->  /ref"
echo "Output dir:     $OUT_DIR   -->  /out"
echo
echo "[COMMAND INSIDE CONTAINER]"
echo "TRsv call -b /data/$(basename "$CRAM") \\"
echo "          -r /ref/$(basename "$REF") \\"
echo "          -x $PLATFORM --build $BUILD \\"
echo "          -p /out/${SAMPLE}_trsv -n $THREADS"
echo "=============================="
echo

singularity exec \
  --bind "${DATA_DIR}:/data,${REF_DIR}:/ref,${OUT_DIR}:/out" \
  --env MULTALIN="/opt/local/tools/multalin" \
  "$SIF" \
  TRsv call \
    -b "/data/$(basename "$CRAM")" \
    -r "/ref/$(basename "$REF")" \
    -x "$PLATFORM" \
    --build "$BUILD" \
    -p "/out/${SAMPLE}_trsv" \
    -n "$THREADS"
