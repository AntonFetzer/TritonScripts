#!/bin/bash
#SBATCH --job-name=uppsala-uniformity-long
#SBATCH --partition=batch-bdw
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/slurm/uniformity-long-%j.out
#SBATCH --error=logs/slurm/uniformity-long-%j.err

set -euo pipefail
module load scicomp-python-env/2025.2

ANALYSIS_DIR="${ANALYSIS_DIR:-$SLURM_SUBMIT_DIR}"
PYTHON_REPO="${PYTHON_REPO:-/scratch/work/fetzera1/Python}"
export PYTHONPATH="${PYTHON_REPO}:${PYTHONPATH:-}"
cd "$ANALYSIS_DIR"

# One sequential multi-minute job avoids excessive scheduler traffic.
for experiment_index in $(seq 0 12); do
  python3 -m UppsalaAcceleratorLogs.time_spatial_uniformity \
    --manifest metadata/manifest.json \
    --results-dir results/final \
    --experiment-index "$experiment_index" \
    --output-dir results/exploratory/time_spatial_uniformity
done
