#!/bin/bash
#SBATCH --job-name=uppsala-parse
#SBATCH --partition=batch-bdw
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
# Four multi-minute tasks avoid excessive scheduler overhead from short jobs.
#SBATCH --array=0-3%4
#SBATCH --output=slurm/parse-%A_%a.out
#SBATCH --error=slurm/parse-%A_%a.err

set -euo pipefail
module load scicomp-python-env/2025.2

ANALYSIS_DIR="${ANALYSIS_DIR:-$SLURM_SUBMIT_DIR}"
PYTHON_REPO="${PYTHON_REPO:-/scratch/work/fetzera1/Python}"
export PYTHONPATH="${PYTHON_REPO}:${PYTHONPATH:-}"
cd "$ANALYSIS_DIR"

python3 -m UppsalaAcceleratorLogs.uppsala_analysis parse \
  --manifest manifest.json \
  --output-dir production/parsed \
  --task-id "$SLURM_ARRAY_TASK_ID" \
  --resume
