#!/bin/bash
#SBATCH --job-name=trsv_call
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

# --- Required ---
CRAM="${1:?Usage: sbatch run_trsv.sbatch /path/to/sample(.haplotagged).cram}"

# --- Core inputs (edit as needed) ---
REF="/pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta"
PLATFORM="ont"       # ont | hifi | clr
BUILD="38"           # 37 | 38 | T2T
THREADS=4
SIF="/projects/ealiyev@xsede.org/TRsv-v1.1.sif"
NO_HOME=0            # set to 1 to pass --no-home

# --- Optional TRsv 'call' extras (leave empty if not used) ---
CONF_FILE=""         # -c /path/to/conf.yaml
REPEAT_BED=""        # -rep /path/to/repeats.bed
REPEAT_U=""          # -reu /path/to/repeats_u.bed
EXCLUDE_BED=""       # -exc /path/to/exclude.bed
LOWCONF_LIST=""      # -lcs /path/to/lowconf.list
GAP_BED=""           # -gb /path/to/gaps.bed
TE_FASTA=""          # -tf /path/to/te.fasta

# --- Derived paths ---
SAMPLE=$(basename "$CRAM" .haplotagged.cram)
SAMPLE="${SAMPLE%.cram}"
OUT_PREFIX="/pl/active/dashnowlab/projects/TR-benchmarking/trsv/${SAMPLE}_trsv"

DATA_DIR=$(dirname "$CRAM")
REF_DIR=$(dirname "$REF")
OUT_DIR=$(dirname "$OUT_PREFIX")

# Per-sample temp dir (like the Perl wrapper)
TMP_DIR="${OUT_DIR}/$(basename "${OUT_PREFIX}").temp"
mkdir -p "$TMP_DIR"

# --- Helper to conditionally add binds ---
add_bind() {
  local p="$1"
  if [[ -n "$p" ]]; then
    local d
    d="$(dirname "$p")"
    if [[ -d "$d" ]]; then
      BIND_LIST+=("$d")
    fi
  fi
}

# --- Build bind list from all needed dirs ---
declare -a BIND_LIST
BIND_LIST+=("$DATA_DIR" "$REF_DIR" "$OUT_DIR" "$TMP_DIR")

add_bind "$CONF_FILE"
add_bind "$REPEAT_BED"
add_bind "$REPEAT_U"
add_bind "$EXCLUDE_BED"
add_bind "$LOWCONF_LIST"
add_bind "$GAP_BED"
add_bind "$TE_FASTA"

# If your host actually has ML training dirs like the Perl wrapper expects,
# bind them; otherwise this silently skips.
if [[ "$PLATFORM" == "ont" && -d "/opt/local/tools/TRsv/Data/ML_train/ONT" ]]; then
  BIND_LIST+=("/opt/local/tools/TRsv/Data/ML_train/ONT")
elif [[ "$PLATFORM" == "clr" && -d "/opt/local/tools/TRsv/Data/ML_train/CLR" ]]; then
  BIND_LIST+=("/opt/local/tools/TRsv/Data/ML_train/CLR")
fi

# De-duplicate bind entries
if command -v python3 >/dev/null 2>&1; then
  # quick uniq while preserving order
  mapfile -t BIND_LIST < <(printf '%s\n' "${BIND_LIST[@]}" | awk '!seen[$0]++')
fi

# Compose --bind argument: host -> inside-container mount points
# We’ll keep your stable mounts and add the rest at same-path mounts.
BIND_SPEC="${DATA_DIR}:/data,${REF_DIR}:/ref,${OUT_DIR}:/out,${TMP_DIR}:/work"
for d in "${BIND_LIST[@]}"; do
  # Skip ones we already mapped explicitly
  [[ "$d" == "$DATA_DIR" || "$d" == "$REF_DIR" || "$d" == "$OUT_DIR" || "$d" == "$TMP_DIR" ]] && continue
  BIND_SPEC+=",${d}"
done

# --- Load Singularity module (best-effort) ---
module load singularity 2>/dev/null || true

# --- Compose TRsv args (NOTE: -p is now relative) ---
TRSV_ARGS=( call
  -b "/data/$(basename "$CRAM")"
  -r "/ref/$(basename "$REF")"
  -x "$PLATFORM"
  --build "$BUILD"
  -p "$(basename "$OUT_PREFIX")"   # e.g., OMLR22_002_trsv (no /out/)
  -n "$THREADS"
)
[[ -n "$CONF_FILE"    ]] && TRSV_ARGS+=( -c  "$(basename "$CONF_FILE")" )
[[ -n "$REPEAT_BED"   ]] && TRSV_ARGS+=( -rep "$(basename "$REPEAT_BED")" )
[[ -n "$REPEAT_U"     ]] && TRSV_ARGS+=( -reu "$(basename "$REPEAT_U")" )
[[ -n "$EXCLUDE_BED"  ]] && TRSV_ARGS+=( -exc "$(basename "$EXCLUDE_BED")" )
[[ -n "$LOWCONF_LIST" ]] && TRSV_ARGS+=( -lcs "$(basename "$LOWCONF_LIST")" )
[[ -n "$GAP_BED"      ]] && TRSV_ARGS+=( -gb  "$(basename "$GAP_BED")" )
[[ -n "$TE_FASTA"     ]] && TRSV_ARGS+=( -tf  "$(basename "$TE_FASTA")" )

# --- Run (NOTE: --pwd /out, and -p is relative) ---
set -x
if [[ "$NO_HOME" -eq 1 ]]; then
  singularity exec --no-home --pwd /out --bind "$BIND_SPEC" --env MULTALIN="/opt/local/tools/multalin" "$SIF" TRsv "${TRSV_ARGS[@]}"
else
  singularity exec --pwd /out --bind "$BIND_SPEC" --env MULTALIN="/opt/local/tools/multalin" "$SIF" TRsv "${TRSV_ARGS[@]}"
fi
set +x

if [[ "$NO_HOME" -eq 1 ]]; then
  singularity exec \
    --no-home \
    --bind "$BIND_SPEC" \
    --env MULTALIN="/opt/local/tools/multalin" \
    "$SIF" \
    TRsv "${TRSV_ARGS[@]}"
else
  singularity exec \
    --bind "$BIND_SPEC" \
    --env MULTALIN="/opt/local/tools/multalin" \
    "$SIF" \
    TRsv "${TRSV_ARGS[@]}"
fi
set +x
