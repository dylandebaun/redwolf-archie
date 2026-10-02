#!/bin/bash
# submit_archie_real.sh — run on the login node from the repository root.
# Submits stages 0-5 as one serial job, then stage 6 as a per-chromosome
# array that starts only if 0-5 succeed, then stage 7 (plots) after the array.
# Stage 8 (apply the trained model) is submitted if a trained model already
# exists in model_dir, or if MODEL_JOB=<train job id> is given (it then waits
# for that job too; see scripts/sim/submit_model.sh).
#
#   bash scripts/submit_archie_real.sh config/rw_hist_safe_53g.yaml
#
# Extra sbatch options (partition, mail, account) can be passed via SBATCH_EXTRA,
# e.g. SBATCH_EXTRA="--mail-type=ALL --mail-user=you@example.org".
set -euo pipefail
CONFIG="${1:?usage: $0 config/<experiment>.yaml}"
cd "$(dirname "${BASH_SOURCE[0]}")/.."
mkdir -p logs

PY="${PYTHON:-python3}"
NCHR=$("${PY}" -c "import yaml,sys;print(len(yaml.safe_load(open(sys.argv[1]))['chromosomes']))" "${CONFIG}")
EXPORT="ALL,CONFIG=${CONFIG}"
[[ -n "${SIM_META:-}" ]] && EXPORT="${EXPORT},SIM_META=${SIM_META}"

J1=$(sbatch --parsable ${SBATCH_EXTRA:-} --job-name=archie_prep \
  --export="${EXPORT},STAGES=0 1 2 3 4 5" scripts/run_archie_real.sh)
J2=$(sbatch --parsable ${SBATCH_EXTRA:-} --job-name=archie_stats \
  --dependency=afterok:${J1} --array=1-${NCHR} --cpus-per-task=1 \
  --export="${EXPORT},STAGES=6" scripts/run_archie_real.sh)
J3=$(sbatch --parsable ${SBATCH_EXTRA:-} --job-name=archie_plot \
  --dependency=afterok:${J2} --cpus-per-task=1 --time=4:00:00 \
  --export="${EXPORT},STAGES=7" scripts/run_archie_real.sh)

MODEL_DIR=$("${PY}" scripts/emit_real_env.py --config "${CONFIG}" | sed -n "s/^MODEL_DIR='\(.*\)'$/\1/p")
if [[ -n "${MODEL_JOB:-}" || -f "${MODEL_DIR}/trained_model.Rdata" ]]; then
  J4=$(sbatch --parsable ${SBATCH_EXTRA:-} --job-name=archie_apply \
    --dependency=afterok:${J2}${MODEL_JOB:+:${MODEL_JOB}} --cpus-per-task=1 --time=4:00:00 \
    --export="${EXPORT},STAGES=8" scripts/run_archie_real.sh)
fi

echo "stages 0-5 : ${J1}"
echo "stage 6    : ${J2} (array 1-${NCHR}, after ${J1})"
echo "stage 7    : ${J3} (after ${J2})"
if [[ -n "${J4:-}" ]]; then
  echo "stage 8    : ${J4} (after ${J2}${MODEL_JOB:+ and model job ${MODEL_JOB}})"
else
  echo "stage 8    : not submitted (no model in ${MODEL_DIR}; run scripts/sim/submit_model.sh or set MODEL_JOB)"
fi
