#!/bin/bash
#SBATCH --job-name=uppsala-uniformity
#SBATCH --partition=batch-bdw
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=0-12%13
#SBATCH --output=slurm/uniformity-%A_%a.out
#SBATCH --error=slurm/uniformity-%A_%a.err

set -euo pipefail
module load scicomp-python-env/2025.2

ANALYSIS_DIR="${ANALYSIS_DIR:-$SLURM_SUBMIT_DIR}"
PYTHON_REPO="${PYTHON_REPO:-/scratch/work/fetzera1/Python}"
export PYTHONPATH="${PYTHON_REPO}:${PYTHONPATH:-}"
cd "$ANALYSIS_DIR"

python3 -m UppsalaAcceleratorLogs.time_spatial_uniformity \
  --manifest manifest.json \
  --results-dir production/results \
  --experiment-index "$SLURM_ARRAY_TASK_ID" \
  --output-dir production/results/time_spatial_uniformity
