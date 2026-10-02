#!/bin/bash
#SBATCH --job-name=archie_sim_train
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --time=60:00:00
#SBATCH --output=logs/slurm-%j_%x.out

# create_training.sh — thin wrapper; all logic lives in create_sims.sh.
#   sbatch --export=ALL,CONFIG=config/x.yaml scripts/sim/create_training.sh   (from the repo root)
set -e
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"
REPO_DIR="${REPO_DIR}" bash "${REPO_DIR}/scripts/sim/create_sims.sh" training_data.txt TRAIN
