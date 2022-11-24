#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=1G
##SBATCH --partition=debug
#SBATCH --output=out.log
#SBATCH --exclusive


srun python ParalelTest.py


## 2000000000
## 100000000
