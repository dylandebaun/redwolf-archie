#!/bin/bash
# train.sh — fit the ArchIE logistic model on a simulation config's training
# data and calibrate thresholds on its test data, using the shared
# scripts/sim/train.R and calibrate_thresholds.R. Run from the repo root:
#
#   bash simulations/train.sh simulations/config/final_53g.yaml
#
# Needs training_data.txt and test_data.txt in the config's out_dir
# (simulate_archie.py run --set train / --set test). R pkgs: data.table, jsonlite.
# Writes <out_dir>/trained_model.Rdata and <out_dir>/thresholds.*.tsv, the
# layout stage 8 of scripts/run_archie_real.sh expects in model_dir.
set -euo pipefail
CONFIG="${1:?usage: $0 simulations/config/<name>.yaml}"
cd "$(dirname "${BASH_SOURCE[0]}")/.."

eval "$(python3 simulations/simulate_archie.py meta --config "${CONFIG}")"
export SIM_META="${MODEL_DIR}/sim_meta.json"

Rscript scripts/sim/train.R

Rscript scripts/sim/calibrate_thresholds.R \
  --train "${MODEL_DIR}/training_data.txt" \
  --test  "${MODEL_DIR}/test_data.txt" \
  --model-rdata "${MODEL_DIR}/trained_model.Rdata" \
  --sim-meta "${SIM_META}" \
  --target-fpr "${TARGET_FPR}" \
  --include-inbetween "${INCLUDE_INBETWEEN}" \
  --out-prefix "${MODEL_DIR}/thresholds"
echo "model: ${MODEL_DIR}/trained_model.Rdata  thresholds: ${MODEL_DIR}/thresholds.*"
