#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash header_fixer.sh /path/to/ont_bam.list /path/to/outdir

LIST="${1:?Provide list file}"
OUTDIR="${2:?Provide output dir}"
mkdir -p "${OUTDIR}"

OLD_REF='/humvar-work-sup/7e/e65de59792c831c956ad96d5d902f0/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta'
NEW_REF='/pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta'

# count non-empty lines
N=$(grep -cv '^\s*$' "${LIST}" || true)
[[ "${N}" -gt 0 ]] || { echo "ERROR: empty list: ${LIST}" >&2; exit 1; }

sbatch \
  --array=1-"${N}" \
  --output="${OUTDIR}/slurm.%x.%A_%a.out" \
  --export=ALL,LIST="${LIST}",OUTDIR="${OUTDIR}",OLD_REF="${OLD_REF}",NEW_REF="${NEW_REF}" \
  urfix_reheader.sbatch
