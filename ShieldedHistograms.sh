#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=75G
#SBATCH --output=ShieldedHistograms.log

srun python ShieldedHistograms.py ProtonsFull 75
srun python ShieldedHistograms.py Protons10MeV 75
srun python ShieldedHistograms.py ElectronsFull 5
srun python ShieldedHistograms.py Electrons500keV 5
