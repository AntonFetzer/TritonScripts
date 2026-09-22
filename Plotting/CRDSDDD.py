"""Aggregate the LED displacement-damage tallies of one CRDS run, any campaign.

The DDD sibling of Plotting/CRDSTID.py, and it reads the SAME result files.
CRDS scores both instruments from one set of primaries, so a run directory
holds TID and DDD blocks in the same CSVs; CRDSTID.py reports the RadFET side
and this reports the LED side. Campaigns are defined in Dependencies/CRDSCampaigns.py.

    python3 Plotting/CRDSDDD.py --campaign HUS --expected-files 100
    python3 Plotting/CRDSDDD.py --campaign Kumpula --expected-files 100 \
        --folder 10MeVProton-CRDS-100umEpoxy
    python3 Plotting/CRDSDDD.py --campaign HUS --expected-files 20 \
        --job-id 20392554 --folder PDD-Electron-CRDS-100mmField-pilot

Displacement damage is quoted in MeV/g, as GRAS writes it. The coefficient is
per incident particle/cm2 because the run normalises in FLUENCE/CURRENT mode,
so a delivered DDD is the coefficient times the delivered fluence. Kumpula has
no delivered fluence in the model yet; pass ``--reference-fluence`` to quote one.

The two tallies are the 350 nm active layer and the substrate under it. DDD is
intensive, MeV per gram, so a field uniform through the die gives the same
answer in both despite the 500x thickness difference. They are reported side by
side for exactly that reason: at HUS they agree, and at Kumpula, where the
protons stop in or near the die, their ratio is what the epoxy sweep tracks.
"""

import argparse
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.CRDSCampaigns import CAMPAIGNS  # noqa: E402
from Dependencies.CRDSDDDReport import reportDDD  # noqa: E402


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

    reportDDD(arguments, parser)


if __name__ == "__main__":
    main()
