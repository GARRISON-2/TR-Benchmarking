#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   build_catalogs.sh "TR_catalog_for_vamos.shard_*_of_05.*.tsv.gz" OUTDIR
#   build_catalogs.sh "TR_catalog_for_vamos.shard_*_of_05.*.tsv.gz" OUTDIR --include bed1.bed bed2.bed ...

SHARD_GLOB="${1:?Provide shard glob (quote it)}"
OUTDIR="${2:?Provide output dir}"
shift 2

INCLUDE_BEDS=()
if [[ "${1:-}" == "--include" ]]; then
  shift
  # collect all remaining args as BEDs
  while [[ $# -gt 0 ]]; do
    INCLUDE_BEDS+=("$1")
    shift
  done
fi

mkdir -p "${OUTDIR}"

module load htslib >/dev/null 2>&1 || true
module load bedtools/2.29.1 >/dev/null 2>&1 || true

command -v bgzip >/dev/null 2>&1 || { echo "ERROR: bgzip not found (htslib)"; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "ERROR: tabix not found (htslib)"; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found"; exit 1; }

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

mapfile -t SHARDS < <(ls -1 ${SHARD_GLOB} 2>/dev/null | sort -V)
[[ "${#SHARDS[@]}" -gt 0 ]] || { echo "ERROR: no files matched: ${SHARD_GLOB}"; exit 1; }
echo "Found ${#SHARDS[@]} shard(s)."

HEADER="$(zcat "${SHARDS[0]}" | head -n 1 || true)"

MERGED_NOHDR="${TMPDIR}/vamos.merged.nohdr.tsv"
SORTED_NOHDR="${TMPDIR}/vamos.sorted.nohdr.tsv"
MASTER_NOHDR="${TMPDIR}/vamos.master.nohdr.tsv"   # sorted, and optionally filtered

# 1) merge (no header/comments)
{
  for f in "${SHARDS[@]}"; do
    zcat "$f" | awk 'NR==1{next} /^#/{next} {print}'
  done
} > "${MERGED_NOHDR}"

# 2) true genomic sort ONCE (chr1..22, X, Y, M/MT, then others)
awk -v OFS="\t" '
function chr_rank(chr,   c, n) {
  c = chr
  sub(/^chr/,"",c)
  if (c ~ /^[0-9]+$/) { n=c+0; if (n>=1 && n<=22) return n }
  if (c=="X")  return 23
  if (c=="Y")  return 24
  if (c=="M" || c=="MT") return 25
  return 99
}
{
  r = chr_rank($1)
  c = $1; sub(/^chr/,"",c)
  cn = (c ~ /^[0-9]+$/) ? c+0 : 999
  print r, cn, $0
}' "${MERGED_NOHDR}" \
| sort -t $'\t' -k1,1n -k2,2n -k3,3 -k4,4n -k5,5n \
| cut -f3- > "${SORTED_NOHDR}"

# 3) OPTIONAL: filter by overlap with one or more BEDs (union)
# bedtools intersect can take multiple -b files: -b bed1 bed2 bed3 ...
if [[ "${#INCLUDE_BEDS[@]}" -gt 0 ]]; then
  echo "Filtering by overlap with ${#INCLUDE_BEDS[@]} BED(s)..."
  # keep records that overlap ANY bed (-u); keep full record (-wa)
  bedtools intersect -wa -u \
    -a "${SORTED_NOHDR}" \
    -b "${INCLUDE_BEDS[@]}" \
    > "${MASTER_NOHDR}"
else
  MASTER_NOHDR="${SORTED_NOHDR}"
fi

# -----------------------------
# Build catalogs from MASTER_NOHDR (already sorted; maybe filtered)
# -----------------------------

# a) Vamos (as-is; can keep extra columns)
VAMOS_BED="${OUTDIR}/catalog.vamos.bed"
cp "${MASTER_NOHDR}" "${VAMOS_BED}"
echo "Wrote: ${VAMOS_BED}"

# b) Strkit: chr start end motif
STRKIT_BED="${OUTDIR}/catalog.strkit.bed"
awk -v OFS="\t" '{print $1,$2,$3,$4}' "${MASTER_NOHDR}" > "${STRKIT_BED}"
echo "Wrote: ${STRKIT_BED}"

# c/d/e) Strdust / Straglr / Medaka: chr start end ID (chromNo-start-end-motif)
STRDUST_BED="${OUTDIR}/catalog.strdust.bed"
awk -v OFS="\t" '
function stripchr(c){ if (c ~ /^chr/) sub(/^chr/,"",c); return c }
{
  c=stripchr($1)
  id=c "-" $2 "-" $3 "-" $4
  print $1,$2,$3,id
}' "${MASTER_NOHDR}" > "${STRDUST_BED}"
echo "Wrote: ${STRDUST_BED}"

STRAGLR_BED="${OUTDIR}/catalog.straglr.bed"
cp "${STRDUST_BED}" "${STRAGLR_BED}"
echo "Wrote: ${STRAGLR_BED}"

MEDAKA_BED="${OUTDIR}/catalog.medaka.bed"
cp "${STRDUST_BED}" "${MEDAKA_BED}"
echo "Wrote: ${MEDAKA_BED}"

# f) LongTR: chr start(1-based) end motif ID (ID uses 0-based coords)
LONGTR_BED="${OUTDIR}/catalog.longtr.bed"
awk -v OFS="\t" '
function stripchr(c){ if (c ~ /^chr/) sub(/^chr/,"",c); return c }
{
  c=stripchr($1)
  id=c "-" $2 "-" $3 "-" $4
  print $1, ($2+1), $3, $4, id
}' "${MASTER_NOHDR}" > "${LONGTR_BED}"
echo "Wrote: ${LONGTR_BED}"

# g) Atarva: chr start end motif motif_size (motif_size = col7)
ATARVA_BED="${OUTDIR}/catalog.atarva.bed"
ATARVA_GZ="${ATARVA_BED}.gz"
awk -v OFS="\t" '{print $1,$2,$3,$4,$7}' "${MASTER_NOHDR}" > "${ATARVA_BED}"
bgzip -f "${ATARVA_BED}"
tabix -f -p bed "${ATARVA_GZ}"
echo "Wrote: ${ATARVA_GZ} and ${ATARVA_GZ}.tbi"

echo "All catalogs written to: ${OUTDIR}"
