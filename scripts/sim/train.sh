#!/bin/bash
#SBATCH --job-name=archie_train
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --time=60:00:00
#SBATCH --output=logs/slurm-%j_%x.out
#
# train.sh — fit the ArchIE logistic regression on training_data.txt, then
# calibrate probability thresholds on test_data.txt for the target FPRs.
#   sbatch --export=ALL,CONFIG=config/x.yaml scripts/sim/train.sh   (from the repo root)
# Writes <model_dir>/trained_model.Rdata and <model_dir>/thresholds.*.tsv
set -eo pipefail
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
cd "${REPO_DIR}"
CONFIG="${CONFIG:?set CONFIG=config/<experiment>.yaml}"
PYTHON="${PYTHON:-python3}"

eval "$("${PYTHON}" scripts/sim/emit_sim_meta.py --config "${CONFIG}")"
export SIM_META

if [[ -n "${CONDA_BASE:-}" ]]; then export PATH="${CONDA_BASE}/bin:${PATH}"; fi
eval "$(conda shell.bash hook)"
conda activate "${ENV_ARCHIE}"

Rscript scripts/sim/train.R

Rscript scripts/sim/calibrate_thresholds.R \
  --train "${MODEL_DIR}/training_data.txt" \
  --test  "${MODEL_DIR}/test_data.txt" \
  --model-rdata "${MODEL_DIR}/trained_model.Rdata" \
  --sim-meta "${SIM_META}" \
  --target-fpr "${TARGET_FPR}" \
  --include-inbetween "${INCLUDE_INBETWEEN}" \
  --out-prefix "${MODEL_DIR}/thresholds"
