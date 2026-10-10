#!/bin/bash
# submit.sh — admixture-fraction sensitivity: for every row of grid.tsv (panel x
# fraction) submit
#   sims (train || test)  ->  train.sh (fit + calibrate)  ->  apply_model.R on the
#   panel's real windows  ;  check_fit.py after the training sims
# then one summary job after everything. Run on the login node from the
# repository root, after make_configs.py:
#
#   python3 sensitivity/admixture/make_configs.py
#   bash sensitivity/admixture/submit.sh
#
# Options (environment): ONLY=<panel> submits one panel only (no summary job);
# SIM_PARTITION=<partition> sends the simulation jobs there (e.g. gpu);
# REAL_DEP=<jobid> makes the apply/fit jobs also wait for that job (e.g. the
# panel's stage-6 feature array).
#
# Outputs per model: sensitivity/admixture/results/<panel>_f<frac>/
#   training_data.txt, test_data.txt, trained_model.Rdata, thresholds.*.tsv,
#   fit/ (check_fit.py), calls/<tag>.{windows.tsv.gz,per_sample.tsv}
# Only files inside sensitivity/admixture/ are written.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
H=sensitivity/admixture
BASE="${CONDA_BASE:-/home/ddebaun/mendel-nas1/miniconda3}"
PY="${BASE}/bin/python3"                      # numpy, scipy, sklearn, pyyaml, pandas, matplotlib
RBIN="${BASE}/envs/archie/bin"               # Rscript + data.table
R=$PWD
mkdir -p "$H/logs"
ALL=()

while IFS=$'\t' read -r panel frac cfg real model_dir; do
  [[ "$panel" == panel ]] && continue
  [[ -n "${ONLY:-}" && "$panel" != "$ONLY" ]] && continue
  SP=${SIM_PARTITION:+-p $SIM_PARTITION}
  RD=${REAL_DEP:+:$REAL_DEP}
  tag=$(basename "$cfg" .yaml)
  # sizes from the first runs: sims peak < 1 GB; historical training ~10 min on 12 CPUs,
  # ~30 min on 4. Small requests backfill much sooner on a busy cluster.
  cpus=4
  SIM="cd $R && $PY simulations/simulate_archie.py run --config $cfg --workers $cpus"
  j1=$(sbatch --parsable -J "sa_simtr_$tag" $SP -c $cpus --mem=4G -t 6:00:00 -o "$H/logs/%x_%j.out" --wrap "$SIM --set train")
  j2=$(sbatch --parsable -J "sa_simte_$tag" $SP -c $cpus --mem=4G -t 3:00:00 -o "$H/logs/%x_%j.out" --wrap "$SIM --set test")
  j3=$(sbatch --parsable -J "sa_train_$tag" -d afterok:$j1:$j2 -c 1 --mem=4G -t 3:00:00 -o "$H/logs/%x_%j.out" \
        --wrap "cd $R && export PATH=$BASE/bin:$RBIN:\$PATH && bash simulations/train.sh $cfg")
  j4=$(sbatch --parsable -J "sa_apply_$tag" -d afterok:$j3$RD -c 1 --mem=8G -t 2:00:00 -o "$H/logs/%x_%j.out" \
        --wrap "cd $R && $RBIN/Rscript scripts/apply_model.R --stats-glob '$real/chr*.stats.tsv.gz' \
                --model-rdata $model_dir/trained_model.Rdata --thresholds $model_dir/thresholds.recommended_by_fpr.tsv \
                --sim-meta $model_dir/sim_meta.json --out-prefix $model_dir/calls/$tag")
  j5=$(sbatch --parsable -J "sa_fit_$tag" -d afterok:$j1$RD -c 2 --mem=6G -t 2:00:00 -o "$H/logs/%x_%j.out" \
        --wrap "cd $R && $PY simulations/check_fit.py --config $cfg --real $real")
  ALL+=("$j3" "$j4" "$j5")
  printf "%-22s sims %s/%s -> train %s -> apply %s ; fit %s\n" "$tag" "$j1" "$j2" "$j3" "$j4" "$j5"
done < "$H/grid.tsv"

if [[ -n "${ONLY:-}" ]]; then echo "ONLY=$ONLY: no summary job submitted"; exit 0; fi
deps=$(IFS=:; echo "${ALL[*]}")
js=$(sbatch --parsable -J sa_summary -d "afterany:$deps" -c 1 --mem=8G -t 2:00:00 -o "$H/logs/%x_%j.out" \
      --wrap "cd $R && $PY $H/summarize.py")
echo "summary (after all models, afterany): $js"
