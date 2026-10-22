"""Pool a mono-energy response-kernel sweep into the kernel CSV the fits read.

    python3 -m RadEx.BuildKernelSummary --sweep <kernel dir> --output <csv> [--channel-offset 12]

Every energy directory listed in the sweep's energies.csv is pooled with
Dependencies.TotalDose.totalDose over all CSV files in its Res/ folder: the dose is
the entry-weighted mean of the files and the error sqrt(sum(N_i^2 sigma_i^2))/sum(N_i).
Files of different length, such as the intermediate snapshot a timed-out task
leaves, therefore enter with their own weight, and totalDose's Birge-ratio check
flags files that differ by more than their seed. totalDose returns kRad; the kernel
is written in rad.

The sweep scores 24 volumes: gate oxides in columns 0-11 and silicon dies in
columns 12-23. --channel-offset selects which block becomes channels 0-11.
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

PYTHON_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_ROOT))
from Dependencies.TotalDose import totalDose  # noqa: E402

KRAD_TO_RAD = 1000.0
N_CHANNELS = 12
FIELDNAMES = ["energy_MeV", "channel", "mean_dose_rad", "combined_error_rad",
              "relative_error", "runs", "total_primaries"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sweep", type=Path, required=True, help="kernel directory with energies.csv")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--channel-offset", type=int, default=0, help="0 for gate oxides, 12 for silicon dies")
    args = parser.parse_args()

    rows = sorted(csv.DictReader(open(args.sweep / "energies.csv", newline="")),
                  key=lambda r: float(r["energy_MeV"]))
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(FIELDNAMES)
        for r in rows:
            res = args.sweep / r["directory"] / "Res"
            pooled = totalDose(str(res))
            block = slice(args.channel_offset, args.channel_offset + N_CHANNELS)
            dose = np.asarray(pooled["dose"], float)[block] * KRAD_TO_RAD
            error = np.asarray(pooled["error"], float)[block] * KRAD_TO_RAD
            runs = len(list(res.glob("*.csv")))
            primaries = int(np.asarray(pooled["entries"])[0])
            for c in range(N_CHANNELS):
                relative = error[c] / dose[c] if dose[c] > 0 else 0.0
                writer.writerow([f"{float(r['energy_MeV']):g}", c, f"{dose[c]:.10e}",
                                 f"{error[c]:.10e}", f"{relative:.10e}", runs, primaries])
    print(f"\nWrote {args.output}: {len(rows)} energies")


if __name__ == "__main__":
    main()

