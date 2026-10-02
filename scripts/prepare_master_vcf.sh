#!/bin/bash
# prepare_master_vcf.sh — harmonise the PLINK-exported master VCF with the
# conventions the ArchIE pipeline expects.
#
#   bash scripts/prepare_master_vcf.sh IN.vcf.gz OUT.vcf.gz REFERENCE.fa.gz
#
# The 53-genome VCF was written by PLINK 1.9, which differs from a GATK VCF in
# three ways that matter here:
#   1. Sample IDs are "FID_IID" (e.g. ALcoy_ALcoy)  -> stripped to IID (ALcoy).
#   2. Chromosomes use PLINK dog codes (1..38, 39=X, 40/41 = MT/unplaced)
#      -> renamed to canFam3.1 names (chr01..chr38, chrX); 40/41 dropped.
#      Beagle maps and the reference FASTA both use the chrNN names.
#   3. REF/ALT are PLINK's major/minor alleles, not the reference allele
#      (INFO/PR = "provisional reference"). ArchIE reads allele 1 as derived,
#      so IFS, S* and private-SNP counts depend on this coding. We swap REF/ALT
#      (and recode GTs) back to canFam3.1 with `bcftools +fixref -m flip-all`.
#      The data are forward-strand, so every mismatch should resolve as a pure
#      swap; the verification below aborts if any strand flip or residual
#      mismatch remains.
#
# Requires bcftools (>=1.10, with plugins). Writes alongside OUT:
#   *.sample_map.tsv   old -> new sample IDs
#   *.prep_report.txt  site counts, swap counts, post-fix REF check
set -euo pipefail

IN="${1:?input VCF}"
OUT="${2:?output VCF (.vcf.gz)}"
FASTA="${3:?reference FASTA (bgzipped + faidx-indexed)}"
REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CHR_MAP="${CHR_MAP:-${REPO_DIR}/config/chr_rename.tsv}"

[[ "${OUT}" == *.vcf.gz ]] || { echo "ERROR: OUT must end in .vcf.gz" >&2; exit 1; }
mkdir -p "$(dirname "${OUT}")"
STEM="${OUT%.vcf.gz}"
SAMPLE_MAP="${STEM}.sample_map.tsv"
REPORT="${STEM}.prep_report.txt"
TMP="${STEM}.tmp.vcf.gz"
trap 'rm -f "${TMP}"' EXIT

# ---- 1. sample map: "X_X" -> "X" (PLINK FID_IID with FID == IID) ------------
bcftools query -l "${IN}" | awk -v OFS='\t' '{
  n = $0; L = length($0)
  if (L % 2 == 1) {
    h = (L - 1) / 2
    if (substr($0, h + 1, 1) == "_" && substr($0, 1, h) == substr($0, h + 2)) n = substr($0, 1, h)
  }
  print $0, n
}' > "${SAMPLE_MAP}"
dups=$(cut -f2 "${SAMPLE_MAP}" | sort | uniq -d)
[[ -z "${dups}" ]] || { echo "ERROR: renamed sample IDs collide: ${dups}" >&2; exit 1; }

# ---- 2. keep mapped chromosomes, rename, drop PR flag, swap to reference ----
KEEP=$(cut -f1 "${CHR_MAP}" | paste -sd,)
bcftools view -t "${KEEP}" -Ou "${IN}" \
  | bcftools annotate --rename-chrs "${CHR_MAP}" -x INFO/PR -Ou \
  | bcftools +fixref - -Oz -o "${TMP}" -- -f "${FASTA}" -m flip-all \
  > "${STEM}.fixref.log" 2>&1

# ---- 3. rename samples, index ------------------------------------------------
bcftools reheader -s "${SAMPLE_MAP}" -o "${OUT}" "${TMP}"
bcftools index -f -t "${OUT}"

# ---- 4. verify ---------------------------------------------------------------
fixes=$(bcftools query -f '%INFO/FIXREF\n' "${OUT}" | sort | uniq -c)
check=$(bcftools +fixref "${OUT}" -- -f "${FASTA}" 2>&1 | grep '^NS')
mismatch=$(awk -F'\t' '$2 ~ /^ref mismatch/ {print $3}' <<< "${check}")
{
  echo "input:  ${IN}"
  echo "output: ${OUT}"
  echo "fasta:  ${FASTA}"
  echo "date:   $(date)"
  echo "bcftools: $(bcftools --version | head -1)"
  echo
  echo "## sites per chromosome (output)"
  bcftools index -s "${OUT}"
  echo
  echo "## INFO/FIXREF actions (none = already matched; swap,GT = REF/ALT swapped + GTs recoded)"
  echo "${fixes}"
  echo
  echo "## post-fix REF check against FASTA"
  echo "${check}"
} > "${REPORT}"

if grep -q 'flip' <<< "${fixes}"; then
  echo "ERROR: fixref complemented the strand at some sites; inspect ${REPORT}" >&2; exit 1
fi
if [[ "${mismatch}" != "0" ]]; then
  echo "ERROR: ${mismatch} sites still mismatch the reference; inspect ${REPORT}" >&2; exit 1
fi
echo "OK: ${OUT} (all REF alleles match ${FASTA}); report: ${REPORT}"
