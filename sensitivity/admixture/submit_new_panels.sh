#!/bin/bash
# Admixture-fraction sensitivity for two more panels, run from redwolf-extra/:
#   hist1975  13 historical (1975) wolves, SAFE excluded  (config/rw_hist1975_53g.yaml)
#   safe2     Cru301 + Cru792 only                        (config/rw_safe2_53g.yaml, from the repo)
# (historical + SAFE together is the existing "historical" panel.)
#
# Both panels phase the same samples as rw_hist_safe_53g (species.redwolf =
# rw1975safe.txt, coyote_bama + coyote.txt), so its phased per-species VCFs are
# linked in and the real-data track starts at stage 4 (merge) -> 5 -> 6.
# Then, per panel: the 8-fraction grid (submit.sh ONLY=<panel>), the pooled
# model (pool.py -> train -> apply / fit), and one summary + figure job.
#
#   python3 sensitivity/admixture/make_configs.py && bash sensitivity/admixture/submit_new_panels.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
H=sensitivity/admixture
BASE="${CONDA_BASE:-/home/ddebaun/mendel-nas1/miniconda3}"
PY="${BASE}/bin/python3"
RBIN="${BASE}/envs/archie/bin"
R=$PWD
PHASED=results/rw_hist_safe_53g/vcfs
mkdir -p "$H/logs"
export CONDA_BASE="$BASE"

declare -A CFG=([hist1975]=config/rw_hist1975_53g.yaml [safe2]=config/rw_safe2_53g.yaml)
TRAIN=(); LAST=()
for panel in hist1975 safe2; do
  cfg=${CFG[$panel]}
  out=$($PY scripts/emit_real_env.py --config "$cfg" | sed -n "s/^OUT_DIR='\(.*\)'$/\1/p")
  mkdir -p "$out/vcfs"
  for sp in redwolf coyote; do
    for ext in vcf.gz vcf.gz.tbi; do
      ln -sfn "$R/$PHASED/${sp}_phased_ALL.$ext" "$out/vcfs/${sp}_phased_ALL.$ext"
    done
  done
  # real data: merge + eigenstrat (stages 4-5), features per chromosome (stage 6)
  jp=$(sbatch --parsable -J "rp_prep_$panel" -p gpu -c 2 --mem=8G -t 4:00:00 -o "$H/logs/%x_%j.out" \
        --export=ALL,CONFIG=$cfg,STAGES="4 5" scripts/run_archie_real.sh)
  js=$(sbatch --parsable -J "rp_stats_$panel" -p gpu -d afterok:$jp --array=1-39 -c 1 --mem=8G -t 12:00:00 \
        -o "$H/logs/%x_%A_%a.out" --export=ALL,CONFIG=$cfg,STAGES=6 scripts/run_archie_real.sh)
  echo "$panel real data: stages 4-5 $jp -> stage 6 array $js"
  # fraction grid; apply/fit also wait for the features
  g=$(ONLY=$panel SIM_PARTITION=gpu REAL_DEP=$js bash $H/submit.sh)
  echo "$g"
  TRAIN+=($(echo "$g" | sed -n 's/.*-> train \([0-9]*\) .*/\1/p'))
  LAST+=($(echo "$g" | sed -n 's/.*-> apply \([0-9]*\) ; fit \([0-9]*\)$/\1 \2/p'))
  declare "REALJOB_$panel=$js"
done

# pooled "all fractions" model per panel (needs the grid simulations)
deps=$(IFS=:; echo "${TRAIN[*]}")
jpool=$(sbatch --parsable -J sa_pool_new -d afterok:$deps -c 1 --mem=8G -t 2:00:00 -o "$H/logs/%x_%j.out" \
        --wrap "cd $R && $PY $H/pool.py hist1975 safe2")
echo "pool: $jpool"
for panel in hist1975 safe2; do
  cfg=$H/configs/${panel}_pooled.yaml; md=$H/results/${panel}_pooled; tag=${panel}_pooled
  real=results/rw_${panel}_53g/data/window_stats; rj=REALJOB_$panel
  j3=$(sbatch --parsable -J "sa_train_$tag" -d afterok:$jpool -c 1 --mem=4G -t 3:00:00 -o "$H/logs/%x_%j.out" \
        --wrap "cd $R && export PATH=$BASE/bin:$RBIN:\$PATH && bash simulations/train.sh $cfg")
  j4=$(sbatch --parsable -J "sa_apply_$tag" -d afterok:$j3:${!rj} -c 1 --mem=8G -t 2:00:00 -o "$H/logs/%x_%j.out" \
        --wrap "cd $R && $RBIN/Rscript scripts/apply_model.R --stats-glob '$real/chr*.stats.tsv.gz' \
                --model-rdata $md/trained_model.Rdata --thresholds $md/thresholds.recommended_by_fpr.tsv \
                --sim-meta $md/sim_meta.json --out-prefix $md/calls/$tag")
  j5=$(sbatch --parsable -J "sa_fit_$tag" -d afterok:$jpool:${!rj} -c 2 --mem=6G -t 2:00:00 -o "$H/logs/%x_%j.out" \
        --wrap "cd $R && $PY simulations/check_fit.py --config $cfg --real $real")
  LAST+=("$j3" "$j4" "$j5")
  echo "$tag: train $j3 -> apply $j4 ; fit $j5"
done

deps=$(IFS=:; echo "${LAST[*]}")
jsum=$(sbatch --parsable -J sa_summary_new -d afterany:$deps -c 1 --mem=8G -t 2:00:00 -o "$H/logs/%x_%j.out" \
       --wrap "cd $R && $PY $H/summarize.py && $PY $H/plot_pooled_ancestry.py")
echo "summary + pooled-ancestry figure: $jsum"
