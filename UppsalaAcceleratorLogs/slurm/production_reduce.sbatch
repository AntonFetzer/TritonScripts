#!/bin/bash
#SBATCH --job-name=uppsala-reduce
#SBATCH --partition=batch-bdw
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/slurm/reduce-%j.out
#SBATCH --error=logs/slurm/reduce-%j.err

set -euo pipefail
module load scicomp-python-env/2025.2

ANALYSIS_DIR="${ANALYSIS_DIR:-$SLURM_SUBMIT_DIR}"
PYTHON_REPO="${PYTHON_REPO:-/scratch/work/fetzera1/Python}"
export PYTHONPATH="${PYTHON_REPO}:${PYTHONPATH:-}"
cd "$ANALYSIS_DIR"

python3 -m UppsalaAcceleratorLogs.uppsala_analysis reduce \
  --manifest metadata/manifest.json \
  --parsed-dir cache/parsed \
  --output-dir results/final
