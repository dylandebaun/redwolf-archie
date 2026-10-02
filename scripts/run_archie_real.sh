#!/bin/bash
#SBATCH --job-name=archie_real
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=72:00:00
#SBATCH --output=logs/slurm-%j_%x.out
#
# run_archie_real.sh — config-driven ArchIE real-data track.
#
#   0  harmonise master VCF (sample IDs, chr names, REF/ALT -> canFam3.1)
#   1  subsample per-species VCFs (+ phasing-only taxa)
#   2  Beagle phasing per species per chromosome
#   3  concatenate per-species phased chromosomes
#   4  drop phasing-only taxa, merge species
#   5  VCF -> Eigenstrat, split admixed/reference panels, split by chr, ranges
#   6  ArchIE feature extraction (one chr per SLURM array task, or serial)
#   7  plot window stats (optional)
#   8  apply the trained model + calibrated thresholds (needs model_dir; see
#      scripts/sim/submit_model.sh). Not in the default STAGES.
#
# Every stage is idempotent (checks its outputs / checkpoint) so a failed run
# can be resubmitted and resumes at the first incomplete stage.
#
# Run from the repository root (SLURM copies the script, so the repo is located
# via the submit directory):
#   sbatch --export=ALL,CONFIG=config/rw_hist_safe_53g.yaml,STAGES="0 1 2 3 4 5" scripts/run_archie_real.sh
#   sbatch --array=1-39 --export=ALL,CONFIG=...,STAGES=6 scripts/run_archie_real.sh
# or just: bash scripts/submit_archie_real.sh config/rw_hist_safe_53g.yaml
#
# Env vars: CONFIG (required), STAGES, FORCE=1 (ignore checkpoints),
#   SIM_META (sim_meta.json from the training track; enables the n_focal gate),
#   CONDA_BASE (miniconda prefix, if conda is not already on PATH),
#   PYTHON (python3 with PyYAML), CALC (calc_stats_window_data.py),
#   BEAGLE_MEM (default 40g).
set -eo pipefail

REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
cd "${REPO_DIR}"
[[ -f scripts/run_archie_real.sh ]] || { echo "ERROR: ${REPO_DIR} is not the repo root; submit from there" >&2; exit 1; }

CONFIG="${CONFIG:?set CONFIG=config/<experiment>.yaml}"
SIM_META="${SIM_META:-}"
STAGES="${STAGES:-0 1 2 3 4 5 6 7}"
FORCE="${FORCE:-0}"
PYTHON="${PYTHON:-python3}"
CALC="${CALC:-${REPO_DIR}/archie/ArchIE/data/calc_stats_window_data.py}"
BEAGLE_MEM="${BEAGLE_MEM:-40g}"
NTHREADS="${SLURM_CPUS_PER_TASK:-8}"

if [[ -n "${CONDA_BASE:-}" ]]; then export PATH="${CONDA_BASE}/bin:${PATH}"; fi
eval "$(conda shell.bash hook)"

mkdir -p logs

# ---- load config + run the validation gate (aborts here on mismatch) ----
eval "$("${PYTHON}" scripts/emit_real_env.py --config "${CONFIG}" ${SIM_META:+--sim-meta "${SIM_META}"})"

VCFS="${OUT_DIR}/vcfs"
DATA="${OUT_DIR}/data"
CHR_DIR="${DATA}/by_chr"
EIGEN_PREFIX="${DATA}/phased_all"
SPLIT_PREFIX="${DATA}/rw_coy_archie"
RANGES_TSV="${DATA}/rw_coy_archie.chrom_ranges.tsv"
WINDOW_STATS="${DATA}/window_stats"
CKPT="${DATA}/.checkpoints"
mkdir -p "${VCFS}" "${CHR_DIR}" "${WINDOW_STATS}" "${CKPT}"

echo "=== ArchIE real track: ${EXPERIMENT_NAME} (config ${CONFIG}) ==="
echo "    n_focal=${N_FOCAL_HAPS} n_ref=${N_REF_HAPS} window=${WINDOW_SIZE} step=${STEP}"
echo "    species: ${SPECIES_LIST} | merge: ${MERGE_ORDER}"
echo "    stages: ${STAGES} (FORCE=${FORCE})"

run_stage() { [[ " ${STAGES} " == *" $1 "* ]]; }
done_file() { echo "${CKPT}/stage_$1.done"; }
is_done()   { [[ "${FORCE}" != "1" && -f "$(done_file "$1")" ]]; }
mark_done() { touch "$(done_file "$1")"; }

# zero-pad chr (X stays "X")
pad_chr() { if [[ "$1" == "X" ]]; then echo "X"; else printf "%02d" "$1"; fi; }

# ============== STAGE 0: harmonise master VCF ================================
# Shared across experiments: skipped whenever the prepared VCF already exists.
if run_stage 0; then
  if [[ -s "${MASTER_VCF}" && -s "${MASTER_VCF}.tbi" && "${FORCE}" != "1" ]]; then
    echo "[0] prepared master VCF exists: ${MASTER_VCF}"
  else
    echo "[0] Harmonise ${MASTER_VCF_RAW} -> ${MASTER_VCF}"
    conda activate "${ENV_BCFTOOLS}"
    bash scripts/prepare_master_vcf.sh "${MASTER_VCF_RAW}" "${MASTER_VCF}" "${REFERENCE_FASTA}"
  fi
else echo "[0] skip"; fi

# ============================ STAGE 1: subsample ============================
if run_stage 1 && ! is_done 1; then
  echo "[1] Subsample per-species VCFs from master"
  conda activate "${ENV_BCFTOOLS}"
  all_samples=$(bcftools query -l "${MASTER_VCF}")
  for sp in ${SPECIES_LIST}; do
    taxa_var="SP_TAXA__${sp}"; taxa="${!taxa_var}"
    phase_var="SP_PHASE__${sp}"; phase="${!phase_var:-}"
    out="${VCFS}/${sp}.vcf.gz"
    if [[ -s "${out}" ]]; then echo "  ${sp}: exists"; continue; fi
    # union analysis taxa + phasing-only taxa (dedup) for richer phasing context
    samp_list="${VCFS}/${sp}.phasingset.txt"
    cat "${TAXA_DIR}/${taxa}" ${phase:+"${TAXA_DIR}/${phase}"} | awk 'NF && !seen[$0]++' > "${samp_list}"
    missing=$(grep -vxF -f <(echo "${all_samples}") "${samp_list}" || true)
    [[ -z "${missing}" ]] || { echo "ERROR: ${sp} samples not in master VCF: ${missing}" >&2; exit 1; }
    echo "  ${sp}: $(wc -l < "${samp_list}") samples (${taxa}${phase:+ + phasing-only ${phase}})"
    bcftools view -S "${samp_list}" "${MASTER_VCF}" -Oz -o "${out}"
    bcftools index -t "${out}"
  done
  mark_done 1
else echo "[1] skip"; fi

# ============================ STAGE 2: phasing ==============================
if run_stage 2 && ! is_done 2; then
  echo "[2] Beagle phasing per species per chromosome"
  conda activate "${ENV_BEAGLE}"
  for sp in ${SPECIES_LIST}; do
    for chr in ${CHROMS}; do
      CHR=$(pad_chr "${chr}")
      out="${VCFS}/${sp}_phased_chr${CHR}"
      if [[ -s "${out}.vcf.gz" ]]; then echo "  ${sp} chr${CHR}: exists"; continue; fi
      echo "  ${sp} chr${CHR}: phasing"
      beagle -Xmx"${BEAGLE_MEM}" \
        gt="${VCFS}/${sp}.vcf.gz" \
        map="${RECOMB_MAP_DIR}/chr${CHR}.beagle.map" \
        out="${out}" \
        chrom="chr${CHR}" \
        nthreads="${NTHREADS}"
    done
  done
  mark_done 2
else echo "[2] skip"; fi

# ===================== STAGE 3: per-species concat ==========================
if run_stage 3 && ! is_done 3; then
  echo "[3] Concatenate per-species phased chromosomes"
  conda activate "${ENV_BCFTOOLS}"
  for sp in ${SPECIES_LIST}; do
    out="${VCFS}/${sp}_phased_ALL.vcf.gz"
    if [[ -s "${out}" ]]; then echo "  ${sp}: concat exists"; continue; fi
    # explicit order from the config (a glob would sort chrX wherever it falls)
    parts=()
    for chr in ${CHROMS}; do parts+=("${VCFS}/${sp}_phased_chr$(pad_chr "${chr}").vcf.gz"); done
    for f in "${parts[@]}"; do bcftools index -f "$f"; done
    bcftools concat "${parts[@]}" -Oz -o "${out}"
    bcftools index -t "${out}"
  done
  mark_done 3
else echo "[3] skip"; fi

# ===================== STAGE 4: merge species ===============================
PHASED_VCF="${VCFS}/${EXPERIMENT_NAME}_phased_ALL.vcf.gz"
if run_stage 4 && ! is_done 4; then
  echo "[4] Merge species into combined phased VCF (order: ${MERGE_ORDER})"
  conda activate "${ENV_BCFTOOLS}"
  merge_inputs=()
  for sp in ${MERGE_ORDER}; do
    full="${VCFS}/${sp}_phased_ALL.vcf.gz"
    phase_var="SP_PHASE__${sp}"
    if [[ -n "${!phase_var:-}" ]]; then
      # drop phasing-only taxa now that phasing is done
      taxa_var="SP_TAXA__${sp}"
      sub="${VCFS}/${sp}_phased_ALL.analysis.vcf.gz"
      bcftools view -S "${TAXA_DIR}/${!taxa_var}" "${full}" -Oz -o "${sub}"
      bcftools index -f -t "${sub}"
      merge_inputs+=("${sub}")
    else
      merge_inputs+=("${full}")
    fi
  done
  bcftools merge "${merge_inputs[@]}" -Oz -o "${PHASED_VCF}"
  bcftools index -f -t "${PHASED_VCF}"
  mark_done 4
else echo "[4] skip"; fi

# ============== STAGE 5: eigenstrat + panel split + ranges ==================
if run_stage 5 && ! is_done 5; then
  echo "[5] VCF->Eigenstrat, split panels, split by chr, build ranges"
  conda activate "${ENV_ARCHIE}"
  python3 scripts/archie_prep/vcf_to_eigenstrat.py \
    --vcf "${PHASED_VCF}" --out-prefix "${EIGEN_PREFIX}"
  python3 scripts/archie_prep/split_admixed_reference.py \
    --eigen-prefix "${EIGEN_PREFIX}" --out-prefix "${SPLIT_PREFIX}" \
    --admixed-csv "${ADMIXED_CSV}" --reference-csv "${REFERENCE_CSV}"
  python3 scripts/archie_prep/split_panels_by_chromosome.py \
    --split-prefix "${SPLIT_PREFIX}" --chr-dir "${CHR_DIR}"
  python3 scripts/archie_prep/build_chrom_ranges.py \
    --chr-dir "${CHR_DIR}" --out-tsv "${RANGES_TSV}"
  if grep -q 9 "${SPLIT_PREFIX}.admixed.geno" "${SPLIT_PREFIX}.reference.geno"; then
    echo "WARNING: missing genotypes (9) in panel .geno files; ArchIE expects complete phased data" >&2
  fi
  mark_done 5
else echo "[5] skip"; fi

# ============== STAGE 6: feature extraction per chromosome ==================
run_calc_one() {
  local CHR="$1" BEGIN="$2" END="$3"
  local out="${WINDOW_STATS}/${CHR}.stats.tsv.gz"
  if [[ -s "${out}" && "${FORCE}" != "1" ]]; then echo "  ${CHR}: exists -> skip"; return 0; fi
  echo "  ${CHR}: ${BEGIN}-${END}"
  python3 "${CALC}" \
    -s "${CHR_DIR}/${CHR}.snp" \
    -i "${SPLIT_PREFIX}.admixed.ind" \
    -a "${CHR_DIR}/${CHR}.admixed.geno" \
    -r "${CHR_DIR}/${CHR}.reference.geno" \
    -c "${CHR}" -b "${BEGIN}" -e "${END}" \
    -w "${WINDOW_SIZE}" -z "${STEP}" \
    > "${WINDOW_STATS}/${CHR}.stats.tsv" \
    2> "${WINDOW_STATS}/${CHR}.stats.log"
  gzip -f "${WINDOW_STATS}/${CHR}.stats.tsv"
}

if run_stage 6 && ! is_done 6; then
  [[ -f "${CALC}" ]] || { echo "ERROR: ${CALC} not found; run archie/setup_archie.sh" >&2; exit 1; }
  conda activate "${ENV_ARCHIE}"
  if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    # array mode: this task handles exactly one line of RANGES_TSV. No stage
    # checkpoint here; idempotency comes from the per-chr output check.
    line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${RANGES_TSV}")
    if [[ -z "${line}" ]]; then
      echo "[6] array id ${SLURM_ARRAY_TASK_ID}: no line in ${RANGES_TSV}; nothing to do"
      exit 0
    fi
    IFS=$'\t' read -r CHR BEGIN END <<< "${line}"
    echo "[6] array id ${SLURM_ARRAY_TASK_ID} -> ${CHR}"
    run_calc_one "${CHR}" "${BEGIN}" "${END}"
    exit 0
  fi
  echo "[6] ArchIE feature extraction per chromosome (serial)"
  while IFS=$'\t' read -r CHR BEGIN END; do
    [[ -z "${CHR}" ]] && continue
    run_calc_one "${CHR}" "${BEGIN}" "${END}"
  done < "${RANGES_TSV}"
  mark_done 6
else echo "[6] skip"; fi

# ============== STAGE 7: plot window stats ==================================
if run_stage 7 && ! is_done 7; then
  echo "[7] Plot window stats"
  conda activate "${ENV_PLOT}"
  python3 scripts/archie_prep/plot_window_stats.py \
    --input-glob "${WINDOW_STATS}/chr*.stats.tsv.gz" \
    --out-dir "${WINDOW_STATS}/plots" \
    --features s_star,min_dist_ref,private_snps \
    --max-points 50000 --dpi 120
  mark_done 7
else echo "[7] skip"; fi

# ============== STAGE 8: apply trained model ================================
if run_stage 8 && ! is_done 8; then
  echo "[8] Apply trained model from ${MODEL_DIR}"
  THR="${MODEL_DIR}/thresholds.recommended_by_fpr.tsv"
  # also accept models trained by the pre-repo scripts (<name>_thresholds.*)
  [[ -f "${THR}" ]] || THR=$(ls "${MODEL_DIR}"/*thresholds.recommended_by_fpr.tsv 2>/dev/null | head -1)
  for f in "${MODEL_DIR}/trained_model.Rdata" "${THR}" "${MODEL_DIR}/sim_meta.json"; do
    [[ -n "${f}" && -f "${f}" ]] || { echo "ERROR: missing model file in ${MODEL_DIR} (${f:-thresholds})" >&2; exit 1; }
  done
  conda activate "${ENV_ARCHIE}"
  Rscript scripts/apply_model.R \
    --stats-glob "${WINDOW_STATS}/chr*.stats.tsv.gz" \
    --model-rdata "${MODEL_DIR}/trained_model.Rdata" \
    --thresholds "${THR}" \
    --sim-meta "${MODEL_DIR}/sim_meta.json" \
    --out-prefix "${OUT_DIR}/calls/${EXPERIMENT_NAME}"
  mark_done 8
else echo "[8] skip"; fi

echo "=== done. checkpoints in ${CKPT}/ ==="
