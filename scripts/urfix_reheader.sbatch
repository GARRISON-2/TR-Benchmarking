#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --job-name=urfix_reheader

set -euo pipefail

: "${LIST:?Need LIST env var (path to bam/cram list)}"
: "${OUTDIR:?Need OUTDIR env var}"
: "${OLD_REF:?Need OLD_REF env var}"
: "${NEW_REF:?Need NEW_REF env var}"

THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Your list has 2 columns: PATH  SEX
LINE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${LIST}" || true)"
IN="$(awk '{print $1}' <<< "${LINE}")"
SEX="$(awk '{print $2}' <<< "${LINE}")"

if [[ -z "${IN}" ]]; then
  echo "ERROR: no input on line ${SLURM_ARRAY_TASK_ID} of ${LIST}" >&2
  exit 1
fi

base="$(basename "${IN}")"
sample="${base%.cram}"
sample="${sample%.bam}"

mkdir -p "${OUTDIR}"

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

header_fixed="${tmpdir}/header_UR_fixed.sam"
out="${OUTDIR}/${sample}.URfix.cram"

echo "[$(date)] Input : ${IN}"
echo "[$(date)] Sex   : ${SEX:-NA}"
echo "[$(date)] Output: ${out}"

samtools view -H "${IN}" \
  | sed "s|${OLD_REF}|${NEW_REF}|g" \
  > "${header_fixed}"

samtools reheader "${header_fixed}" "${IN}" > "${out}"
samtools index -@ "${THREADS}" "${out}"

echo "[$(date)] Done."
