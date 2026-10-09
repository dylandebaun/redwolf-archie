#!/bin/bash
#SBATCH --job-name=ibdmix
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --array=1-38
#SBATCH --output=logs/slurm-%A_%a_%x.out
#
# run_ibdmix.sbatch — IBDmix with coyote donor(s) ("archaic" role) and 26 red
# wolves as the modern panel, one autosome per array task.
# Adapted from methods_v2/ibdmix_mendel/run_ibdmix.sbatch (Heliconius).
#
#   modern panel = taxa/ibdmix_modern.txt: 13 historical (1975) + 2 SAFE +
#                  11 LSUMZ museum red wolves (9044, dog-like, excluded).
#                  IBDmix estimates modern allele frequencies from all of
#                  them (no -s), and every one gets segments reported.
#   donors       = DONORS (default ALcoy). Each donor genome gets its own
#                  generate_gt + ibdmix run (several donors in one archaic VCF
#                  would be treated as modern samples). Combine afterwards,
#                  e.g. max LOD per site (scripts/ibdmix/combine_donors.py).
#   input        = stage-0 VCF (canFam3.1 REF, unphased; IBDmix uses
#                  unphased genotypes), restricted to the modern panel + all
#                  donors and to biallelic SNPs polymorphic among them (one
#                  site set shared by every donor run). Chromosomes are
#                  renamed chrNN -> N because IBDmix expects integer names.
#   chrX is not run: males are hemizygous there and IBDmix assumes diploids.
#   -d 0.0       report every segment with LOD > 0; apply LOD >= 4 afterwards.
#   -t           add region-level statistics columns.
#   MASK_DIR     optional: folder with <donor>.bed (integer chromosome names);
#                passed as ibdmix -r, which sets LOD to 0 inside those regions
#                (used for each donor's runs of homozygosity).
#
# Run from the repository root:
#   sbatch scripts/ibdmix/run_ibdmix.sbatch                                  # ALcoy
#   sbatch --export=ALL,DONORS="ALcoy ILcoy FLcoy OHcoy",OUT=results/ibdmix_coyotes_53g \
#          scripts/ibdmix/run_ibdmix.sbatch
# Output: <OUT>/<donor>/ibdmix.chrNN.out (one per donor and autosome). The first
# ALcoy-only run (results/ibdmix_alcoy_53g/ibdmix.chrNN.out) predates the
# per-donor subfolders.
set -euo pipefail

REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
cd "${REPO_DIR}"
IBDMIX="${IBDMIX:-/home/ddebaun/mendel-nas1/methods_v2/IBDmix/build/src}"   # generate_gt, ibdmix (commit 6a9d4bc)
VCF="${VCF:-data/RedWolfProject_53genomes_geno10_maf1.canfam31.vcf.gz}"
MODERN="${MODERN:-taxa/ibdmix_modern.txt}"   # e.g. taxa/ibdmix_modern_latx.txt (+ 19 LA/TX hybrids)
DONORS="${DONORS:-ALcoy}"
OUT="${OUT:-results/ibdmix_alcoy_53g}"
MASK_DIR="${MASK_DIR:-}"

if [[ -n "${CONDA_BASE:-}" ]]; then export PATH="${CONDA_BASE}/bin:${PATH}"; fi
eval "$(conda shell.bash hook)"
conda activate genome-assembly       # bcftools

N="${SLURM_ARRAY_TASK_ID:?run as an array job (1-38)}"
CHR=$(printf "chr%02d" "${N}")
WORK="${OUT}/work_${CHR}"
mkdir -p "${OUT}" "${WORK}"

# 1. this chromosome, modern panel + donors, biallelic SNPs polymorphic among them; chrNN -> N
printf '%s\t%s\n' "${CHR}" "${N}" > "${WORK}/rename.txt"
{ cat "${MODERN}"; printf '%s\n' ${DONORS}; } > "${WORK}/samples.txt"
bcftools view -r "${CHR}" -S "${WORK}/samples.txt" -Ou "${VCF}" \
  | bcftools view -m2 -M2 -v snps -c1:minor -Ou \
  | bcftools annotate --rename-chrs "${WORK}/rename.txt" -Oz -o "${WORK}/bi.vcf.gz"
bcftools index -t "${WORK}/bi.vcf.gz"
# SNPs per 1 Mb in this site set (for the map's low-data mask)
bcftools query -f '%POS\n' "${WORK}/bi.vcf.gz" | awk -v C="${CHR}" -v B=1000000 '{n[int($1/B)*B]++} END{for(k in n) print C"\t"k"\t"n[k]}' \
  | sort -k2,2n > "${OUT}/snps_${CHR}.tsv"

# 2. modern VCF (uncompressed: generate_gt needs plain text)
bcftools view -S "${MODERN}" -Ov -o "${WORK}/mod.vcf" "${WORK}/bi.vcf.gz"

# 3. one generate_gt + ibdmix run per donor genome
for DONOR in ${DONORS}; do
  if [[ "${DONORS}" == "ALcoy" && "${OUT}" == "results/ibdmix_alcoy_53g" ]]; then
    RES="${OUT}/ibdmix.${CHR}.out"                   # layout of the original ALcoy-only run
  else
    mkdir -p "${OUT}/${DONOR}"; RES="${OUT}/${DONOR}/ibdmix.${CHR}.out"
  fi
  bcftools view -s "${DONOR}" -Ov -o "${WORK}/arch.vcf" "${WORK}/bi.vcf.gz"
  "${IBDMIX}/generate_gt" -a "${WORK}/arch.vcf" -m "${WORK}/mod.vcf" -o "${WORK}/combined.gt"
  MASK=()
  if [[ -n "${MASK_DIR}" ]]; then
    [[ -f "${MASK_DIR}/${DONOR}.bed" ]] || { echo "ERROR: no mask ${MASK_DIR}/${DONOR}.bed" >&2; exit 1; }
    MASK=(-r "${MASK_DIR}/${DONOR}.bed")
  fi
  "${IBDMIX}/ibdmix" -g "${WORK}/combined.gt" -n "${DONOR}" -o "${RES}" -d 0.0 -t "${MASK[@]}"
  echo "${CHR} ${DONOR}: $(bcftools index -n "${WORK}/bi.vcf.gz") SNPs, $(($(wc -l < "${RES}") - 1)) segments"
done
rm -rf "${WORK}"
echo IBDMIX_DONE
