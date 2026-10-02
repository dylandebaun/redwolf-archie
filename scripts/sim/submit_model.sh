#!/bin/bash
# submit_model.sh — run on the login node from the repository root.
# Renders ms.sh from the config, simulates training + test data in parallel,
# then trains the model and calibrates thresholds once both finish.
#
#   bash scripts/sim/submit_model.sh config/rw_hist_safe_53g.yaml
#
# Prints the training job id; pass it as MODEL_JOB to submit_archie_real.sh so
# the real-data model application (stage 8) waits for it.
set -euo pipefail
CONFIG="${1:?usage: $0 config/<experiment>.yaml}"
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
mkdir -p logs
PY="${PYTHON:-python3}"

[[ -x archie/ArchIE/msmodified/ms ]] || { echo "ERROR: ms not built; run archie/setup_archie.sh" >&2; exit 1; }
eval "$("${PY}" scripts/sim/emit_sim_meta.py --config "${CONFIG}")"
TEMPLATE=$("${PY}" -c "import yaml,sys;print(yaml.safe_load(open(sys.argv[1]))['sim']['demography_template'])" "${CONFIG}")
"${PY}" scripts/sim/render_ms.py --config "${CONFIG}" --template "${TEMPLATE}" --out "${MODEL_DIR}/ms.sh"

EXPORT="ALL,CONFIG=${CONFIG}"
J1=$(sbatch --parsable ${SBATCH_EXTRA:-} --export="${EXPORT}" scripts/sim/create_training.sh)
J2=$(sbatch --parsable ${SBATCH_EXTRA:-} --export="${EXPORT}" scripts/sim/create_test.sh)
J3=$(sbatch --parsable ${SBATCH_EXTRA:-} --dependency=afterok:${J1}:${J2} --export="${EXPORT}" scripts/sim/train.sh)

echo "training sims : ${J1}"
echo "test sims     : ${J2}"
echo "train+calib   : ${J3} (after ${J1}, ${J2})"
echo "model_dir     : ${MODEL_DIR}"
echo
echo "next: MODEL_JOB=${J3} bash scripts/submit_archie_real.sh ${CONFIG}"
