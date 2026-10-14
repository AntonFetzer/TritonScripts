"""Generate the mono-energetic response-kernel sweep on the new RadFET geometry.

Mirrors EnergyKernel-FinalGeometry exactly -- same 60-point energy grid, same 10
replicates of 4,000,000 protons per energy, same mono-energetic source macro --
but points at RadEx-Uppsala-measuredThick-NewRadFET.gdml, confines a 1 um
production cut to the twelve packages, and scores 24 volumes: gate oxides in
columns 0-11 and silicon dies in columns 12-23.

Run from the Python repository root (python3 -m RadExUppsala.GenerateNewRadFETKernel).
Writes nothing outside GRAS/RadEx/RadEx-Uppsala/EnergyKernel-NewRadFET.
"""

import csv
from pathlib import Path
import shutil
import sys

CAMPAIGN = Path("/scratch/work/fetzera1/GRAS/RadEx/RadEx-Uppsala")
OLD_SWEEP = CAMPAIGN / "EnergyKernel-FinalGeometry"
NEW_SWEEP = CAMPAIGN / "EnergyKernel-NewRadFET"

PRIMARIES_PER_TASK = 4000000
REPLICATES = 10

REGION_BLOCK = """
############# Production cuts #############
# A 0.7 mm range cut is meaningless across a 400 nm layer. The fine cut is
# confined to the twelve packages, so the degrader, the air path, the plates,
# the inserts and both PCBs keep the Geant4 default.
# Logical volumes exist only after /run/initialize, so this block must follow it.
/gras/physics/region/addRegion RadFET
/control/loop ../../../addRadFETToRegion.mac N 0 11 1

# Two commands, and the order matters. GRAS creates the region without a
# G4ProductionCuts object, so the stock Geant4 command segfaults if called
# first; GRAS setRegionCut cannot be used alone because it rebuilds that object
# on every call, leaving every particle but the last at 0 fm, and rejects
# proton. Seed the object with GRAS, then let Geant4 set all four cuts at once.
/gras/physics/region/setRegionCut RadFET e- {Cut} mm
/run/setCutForRegion RadFET {Cut} mm
/gras/physics/region/listRegions
"""

DETECTORS = """/gras/geometry/util/listPhysicalVolumes

/gras/analysis/dose/addModule doseDetector

# Columns 0-11 are the gate oxides, channel by channel; columns 12-23 are the
# silicon dies behind them. The oxide is the quantity the RadFET calibration
# refers to; the die is the like-for-like comparison against the legacy
# Si_vol_<N> kernel, and carries about five times better statistics.
/control/loop ../../../addVolumeVT01Gox.mac N 0 11 1
/control/loop ../../../addVolumeVT01Die.mac N 0 11 1

/gras/analysis/dose/doseDetector/setUnit rad
/gras/analysis/dose/doseDetector/tallyIndividualVolumes

/gras/analysis/describe 2
"""

RUN_SH = """#!/bin/bash

#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --output=log/%x-job%A-%a.log
#SBATCH --array=1-{replicates}

sleep $((SLURM_ARRAY_TASK_ID*2))

echo "$(date)"
echo "Cluster: $SLURM_CLUSTER_NAME"
echo "Job Name: $SLURM_JOB_NAME"
echo "Job ID: $SLURM_JOBID"
echo "On Node: $SLURMD_NODENAME"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"

srun gras A.mac
"""

SWEEP_SH = """#!/bin/bash

#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --output=sweep-log/%x-job%A-%a.log
#SBATCH --array=0-{last}%100

sleep $((SLURM_ARRAY_TASK_ID*2))

echo "$(date)"
echo "Cluster: $SLURM_CLUSTER_NAME"
echo "Job Name: $SLURM_JOB_NAME"
echo "Job ID: $SLURM_JOBID"
echo "On Node: $SLURMD_NODENAME"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"

energy_index=$((SLURM_ARRAY_TASK_ID / {replicates}))
replicate=$((SLURM_ARRAY_TASK_ID % {replicates} + 1))
directory=$(sed -n "$((energy_index + 1))p" energy-directories.txt)

echo "Energy directory: $directory"
echo "Replicate: $replicate"

cd "$directory"
srun gras A.mac
"""

A_MAC = """# RadEx Uppsala mono-energy response-kernel point: {energy} MeV.
# New RadFET package module with the production cut confined to the packages.
/control/alias GDML ../../RadEx-Uppsala-measuredThick-NewRadFET
/control/alias NumPart {primaries}
/control/alias Cut 0.001
/control/alias BeamEnergy {energy}
/control/alias Spectrum Clinical/Skandion/MonoEnergyKernel
/control/execute ../EnergyKernel-NewRadFET-main.mac
"""


def build_main_macro():
    source = (OLD_SWEEP / "EnergyKernel-main.mac").read_text(encoding="utf-8")

    replacements = [
        ("# Uppsala entry point. This retains the shared RadEx physics, detector and",
         "# Response-kernel entry point on the new RadFET package geometry. Physics,"),
        ("# analysis configuration while selecting the orientation-correct beam macro.",
         "# beam, normalisation and analysis match EnergyKernel-main.mac; the RadFET\n"
         "# representation, the detector list and the production cuts differ."),
        ("/run/verbose 0", "/run/verbose 1"),
        ("/control/execute ../EnergyKernel-detectors.mac",
         "/control/execute ../EnergyKernel-NewRadFET-detectors.mac"),
        ("/run/initialize", "/run/initialize\n" + REGION_BLOCK),
    ]
    text = source
    for old, new in replacements:
        if text.count(old) != 1:
            raise RuntimeError(
                f"anchor {old[:50]!r} appears {text.count(old)}x, expected once"
            )
        text = text.replace(old, new)
    return text


def main() -> None:
    if not OLD_SWEEP.is_dir():
        sys.exit(f"Reference sweep not found: {OLD_SWEEP}")
    if NEW_SWEEP.exists():
        sys.exit(f"Refusing to overwrite an existing {NEW_SWEEP}")

    with (OLD_SWEEP / "energies.csv").open(encoding="utf-8") as handle:
        grid = list(csv.DictReader(handle))
    if not grid:
        sys.exit("energies.csv is empty")

    NEW_SWEEP.mkdir(parents=True)
    (NEW_SWEEP / "sweep-log").mkdir()

    (NEW_SWEEP / "EnergyKernel-NewRadFET-main.mac").write_text(
        build_main_macro(), encoding="utf-8", newline="\n")
    (NEW_SWEEP / "EnergyKernel-NewRadFET-detectors.mac").write_text(
        DETECTORS, encoding="utf-8", newline="\n")

    directories = []
    rows = []
    for point in grid:
        energy = point["energy_MeV"]
        directory = point["directory"]
        directories.append(directory)

        target = NEW_SWEEP / directory
        (target / "log").mkdir(parents=True)
        (target / "Res").mkdir()
        (target / "A.mac").write_text(
            A_MAC.format(energy=energy, primaries=PRIMARIES_PER_TASK),
            encoding="utf-8", newline="\n")
        (target / "run.sh").write_text(
            RUN_SH.format(replicates=REPLICATES), encoding="utf-8", newline="\n")

        rows.append({
            "energy_MeV": energy,
            "directory": directory,
            "array_tasks": REPLICATES,
            "primaries_per_task": PRIMARIES_PER_TASK,
            "total_primaries": REPLICATES * PRIMARIES_PER_TASK,
        })

    (NEW_SWEEP / "energy-directories.txt").write_text(
        "\n".join(directories) + "\n", encoding="utf-8", newline="\n")

    with (NEW_SWEEP / "energies.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    (NEW_SWEEP / "run-energy-kernel-sweep.sh").write_text(
        SWEEP_SH.format(last=len(grid) * REPLICATES - 1, replicates=REPLICATES),
        encoding="utf-8", newline="\n")

    print(f"Energies: {len(grid)}")
    print(f"Tasks: {len(grid) * REPLICATES}")
    print(f"Primaries per energy: {REPLICATES * PRIMARIES_PER_TASK:,}")
    print(f"Total primaries: {len(grid) * REPLICATES * PRIMARIES_PER_TASK:,}")
    print(f"Written to {NEW_SWEEP}")


if __name__ == "__main__":
    main()
