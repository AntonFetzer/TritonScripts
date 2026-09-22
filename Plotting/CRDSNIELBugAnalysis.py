"""Compatibility entry point for the CRDS NIEL transport diagnostic."""

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Tests.CRDSNIELBugAnalysis import validate  # noqa: E402

if __name__ == "__main__":
    validate()
