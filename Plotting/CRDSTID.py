"""Aggregate the RadFET TID tallies of one CRDS production run, any campaign.

One VT01 RadFET package on the 1.6 mm carrier PCB. Two tallies, the 400 nm
gate oxide and the silicon die behind it, in the order
``CRDS/CRDSRadFETLEDDetector.mac`` fixes. The campaign, its runs and its
reference fluence are defined in Dependencies/CRDSCampaigns.py.

    python3 Plotting/CRDSTID.py --campaign HUS --expected-files 100
    python3 Plotting/CRDSTID.py --campaign Kumpula --expected-files 100 \
        --folder 10MeVProton-CRDS-100umEpoxy

``--job-id`` defaults to the job recorded for a listed run and is required for
an unlisted folder such as a pilot.

HUS doses are quoted at 2e12 electrons/cm2, the assumed reference fluence the
RadEx-HUS results use; that is not an MU calibration and not a measured
delivered dose. The 20 x 20 mm HUS field is below lateral scatter equilibrium,
so see Plotting/CRDSSeries.py before quoting it for a broad-field exposure.

Kumpula has no delivered fluence in the model yet, so its results stay per
proton/cm2 unless ``--reference-fluence`` supplies one.
"""

import argparse
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.CRDSCampaigns import CAMPAIGNS  # noqa: E402
from Dependencies.CRDSTIDReport import reportTID  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--campaign", required=True, choices=list(CAMPAIGNS))
    parser.add_argument("--expected-files", type=int, required=True)
    parser.add_argument("--folder",
                        help="run directory; defaults to the campaign's "
                             "production run")
    parser.add_argument("--job-id", type=int)
    parser.add_argument("--reference-fluence", type=float,
                        help="particles/cm2; defaults to the campaign's")
    arguments = parser.parse_args()

    reportTID(arguments, parser)


if __name__ == "__main__":
    main()
