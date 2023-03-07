#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=1G
##SBATCH --partition=debug
#SBATCH --output=log/Res.log

echo $(date)
echo "Cluster: " $SLURM_CLUSTER_NAME
echo "Job Name: " $SLURM_JOB_NAME
echo "Job ID: " $SLURM_JOBID
echo "On Node: " $SLURMD_NODENAME

srun python3 -u LingyunSim2.py
