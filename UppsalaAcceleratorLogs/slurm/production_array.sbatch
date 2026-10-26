#!/bin/bash
#SBATCH --job-name=uppsala-parse
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
# Four multi-minute tasks avoid excessive scheduler overhead from short jobs.
#SBATCH --array=0-3%4
#SBATCH --output=logs/slurm/parse-%A_%a.out
#SBATCH --error=logs/slurm/parse-%A_%a.err

set -euo pipefail
module load scicomp-python-env/2025.2

ANALYSIS_DIR="${ANALYSIS_DIR:-$SLURM_SUBMIT_DIR}"
PYTHON_REPO="${PYTHON_REPO:-/scratch/work/fetzera1/Python}"
export PYTHONPATH="${PYTHON_REPO}:${PYTHONPATH:-}"
cd "$ANALYSIS_DIR"

python3 -m UppsalaAcceleratorLogs.uppsala_analysis parse \
  --manifest metadata/manifest.json \
  --output-dir cache/parsed \
  --task-id "$SLURM_ARRAY_TASK_ID" \
  --resume \
  "$@"
# Extra arguments restrict or override the selection, e.g. re-parse the RadEx
# part records in place (adds heatmap_positions_c to their cache):
#   sbatch production_array.sbatch \
#     --experiments exp_01,exp_02,exp_03,exp_04,exp_05,exp_06,exp_07,exp_08,exp_09 \
#     --kind part --no-resume
