#!/bin/bash
# create_sims.sh — shared core for training/test data generation.
# Usage (from the repository root):
#   CONFIG=config/x.yaml bash scripts/sim/create_sims.sh training_data.txt TRAIN
#   CONFIG=config/x.yaml bash scripts/sim/create_sims.sh test_data.txt     TEST
#
# PREFIX qualifies the out.* work dir and sims staging dir so TEST and TRAIN
# can run concurrently without clobbering each other.
# Reads everything from the config via emit_sim_meta.py; writes into model_dir.
set -euo pipefail

OUT_FILE_NAME="${1:?usage: create_sims.sh <output_file> [PREFIX]}"
SIM_PREFIX="${2:-}"
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
cd "${REPO_DIR}"
[[ -f scripts/sim/create_sims.sh ]] || { echo "ERROR: ${REPO_DIR} is not the repo root; submit from there" >&2; exit 1; }
CONFIG="${CONFIG:?set CONFIG=config/<experiment>.yaml}"
PYTHON="${PYTHON:-python3}"
CALC="${CALC:-${REPO_DIR}/archie/ArchIE/simulations/calc_stats_ms_rw.py}"

# pull NREPS, WINDOW_SIZE, MODEL_DIR, etc. into the shell + (re)write sim_meta.json
eval "$("${PYTHON}" scripts/sim/emit_sim_meta.py --config "${CONFIG}")"
MS="${MS:-${MODEL_DIR}/ms.sh}"
DIR="${MODEL_DIR}"

if [[ -n "${CONDA_BASE:-}" ]]; then export PATH="${CONDA_BASE}/bin:${PATH}"; fi
eval "$(conda shell.bash hook)"
conda activate "${ENV_ARCHIE}"

[[ -f "${CALC}" ]] || { echo "ERROR: ${CALC} not found; run archie/setup_archie.sh" >&2; exit 1; }
[[ -f "${MS}" ]] || { echo "ERROR: ${MS} not found; run scripts/sim/render_ms.py (submit_model.sh does)" >&2; exit 1; }

# select nreps based on which set we're building
case "${SIM_PREFIX}" in
  TRAIN) NREPS="${NREPS_TRAIN}" ;;
  TEST)  NREPS="${NREPS_TEST}" ;;
  *)     NREPS="${NREPS_TRAIN}" ;;
esac

# per-prefix work dir isolates out.* between concurrent runs
WORK_DIR="${DIR}${SIM_PREFIX:+/${SIM_PREFIX}}"
SIMS_DIR="${WORK_DIR}/sims"
OUT_FILE="${DIR}/${OUT_FILE_NAME}"

mkdir -p "${WORK_DIR}" "${SIMS_DIR}"

echo "Config: nreps=${NREPS} (train=${NREPS_TRAIN} test=${NREPS_TEST}) window=${WINDOW_SIZE} n_focal=${N_FOCAL_HAPS} feature_count=${FEATURE_COUNT} model_dir=${DIR} prefix=${SIM_PREFIX:-<none>}"

# sanity: the rendered ms.sh must match the focal/ref counts in the config
if ! grep -q -- "-I 3 ${N_FOCAL_HAPS} 1 ${N_REF_HAPS}" "${MS}"; then
  echo "ERROR: ${MS} does not contain '-I 3 ${N_FOCAL_HAPS} 1 ${N_REF_HAPS}'." >&2
  echo "       Re-run render_ms.py — ms.sh and config are out of sync." >&2
  exit 1
fi

MS_ABS="$(realpath "${MS}")"

for i in $(seq 1 "${NREPS}"); do
  echo "Simulation ${i}"
  rm -f "${WORK_DIR}"/out.*                 # never reuse a previous replicate's output
  (cd "${WORK_DIR}" && bash "${MS_ABS}")   # writes out.* into WORK_DIR
  python "${CALC}" \
    -s "${WORK_DIR}/out.snp" \
    -a "${WORK_DIR}/out.ADMIXED.geno" \
    -r "${WORK_DIR}/out.2.geno" \
    --anc "${WORK_DIR}/out.ADMIXED.anc" \
    -n "${WINDOW_SIZE}" \
    > "${SIMS_DIR}/${i}.txt"
done

# numeric order so replicate order is stable
for i in $(seq 1 "${NREPS}"); do cat "${SIMS_DIR}/${i}.txt"; done > "${OUT_FILE}"
rm -f "${SIMS_DIR}"/*.txt

# quick column-count check against the contract
NCOL=$(awk '{print NF; exit}' "${OUT_FILE}")
echo "Wrote ${OUT_FILE}: first-row columns=${NCOL} (expected total_cols=${TOTAL_COLS})"
if [[ "${NCOL}" -ne "${TOTAL_COLS}" ]]; then
  echo "WARNING: column count ${NCOL} != expected ${TOTAL_COLS}. Check calc_stats / N." >&2
fi
