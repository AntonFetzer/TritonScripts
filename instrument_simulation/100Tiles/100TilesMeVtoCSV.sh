#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --mem=100G
#SBATCH --output=100TilesMeVtoCSV.log


srun python 100TilesMeVtoCSV.py

