"""Compare the lying and the upright LED on the first Uppsala board.

Two identical GaAs LEDs 10 mm apart in a flat part of the 5.2 x 5.4 cm field,
one lying on its back and one standing upright on its leads. Reports both
LEDs' DDD and substrate TID coefficients, their upright/lying ratios, and the
delivered values from the local fluence over both exposures. The run and its
tally list are UPPSALA_LED_PAIR in Dependencies/CRDSCampaigns.py; the logic is
Dependencies/CRDSLEDPairReport.py.

    python3 Plotting/CRDSLEDPair.py

Writes, beside the simulation data, the same summaries as the other two
Uppsala runs: TotalDose_<folder>.csv and TotalNID_<folder>.csv in the run
folder, and DeliveredDose_CRDS-UppsalaPCB1LEDs.csv in the format of
DeliveredDose_CRDS-Uppsala.csv; plus LEDOrientation_CRDS-Uppsala.csv, the
upright-against-lying comparison with units.
"""

import argparse
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.CRDSCampaigns import UPPSALA_LED_PAIR  # noqa: E402
from Dependencies.CRDSLEDPairReport import reportLEDPair  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.parse_args()
    print()
    for written in reportLEDPair(UPPSALA_LED_PAIR):
        print(f"Wrote {written}")


if __name__ == "__main__":
    main()
