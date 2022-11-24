#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem=0
#SBATCH --partition=debug
#SBATCH --output=Parallel100TilesMeVtoCSV.log
#SBATCH --exclusive

srun python Parallel100TilesMeVtoCSV.py
