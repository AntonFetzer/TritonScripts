"""Per-energy run counts in RadExUppsala.BuildKernelSummary.

A topped-up sweep holds more result files at some energies than at others.
The run count each energy must hold comes from energies.csv unless a fixed
count is forced, and the pooled error and the empirical scatter use that
energy's own count. Result files are built from the committed ReplicaRun
fixture, widened to the 24 volumes the kernel sweeps score, and read back
through Read.ReadDose by the script itself.
"""

import csv
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

from RadExUppsala import BuildKernelSummary as bks

TEMPLATE = Path(__file__).resolve().parent / "Fixtures" / "ReplicaRun" / "TID_replica1.csv"
N_VOLUMES = 24
PRIMARIES = 1000

# Doses and errors in rad per run; Read.ReadDose returns kRad and the script converts back.
RUNS = {
    "001MeVProton": [(1.0, 0.1), (3.0, 0.2)],
    "002MeVProton": [(2.0, 0.1), (4.0, 0.3), (9.0, 0.2)],
}


def write_result(path, dose_rad, error_rad, entries):
    """The template file with N_VOLUMES volumes all carrying the same values."""
    lines = TEMPLATE.read_text(encoding="utf-8").splitlines()
    names = next(i for i, line in enumerate(lines) if "_PV'" in line)
    first_value = next(i for i, line in enumerate(lines) if line.strip()[:1].isdigit())
    n_template = sum(1 for line in lines if line.strip()[:1].isdigit())
    prefix = lines[names].split(",'")[0]
    volume_names = prefix + "," + ",".join(f"'vol_{i}_PV'" for i in range(N_VOLUMES))
    value = f"     {dose_rad:g},       {error_rad:g},       {entries:g},         100"
    out = (lines[:names] + [volume_names] + lines[names + 1:first_value]
           + [value] * N_VOLUMES + lines[first_value + n_template:])
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def make_sweep(root, manifest_runs, with_counts=True):
    rows = []
    for index, (directory, runs) in enumerate(RUNS.items(), start=1):
        res = root / directory / "Res"
        res.mkdir(parents=True)
        for j, (dose, error) in enumerate(runs):
            write_result(res / f"TID_{index}_{j}.csv", dose, error, PRIMARIES)
        row = {"energy_MeV": index, "directory": directory}
        if with_counts:
            row.update(array_tasks=manifest_runs[directory], primaries_per_task=PRIMARIES,
                       total_primaries=manifest_runs[directory] * PRIMARIES)
        rows.append(row)
    with (root / "energies.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class PerEnergyRunCounts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_counts_come_from_the_manifest(self):
        make_sweep(self.root, {"001MeVProton": 2, "002MeVProton": 3})
        energies, means, errors, sems, runs, primaries = bks.aggregate(self.root, 0)

        np.testing.assert_array_equal(energies, [1.0, 2.0])
        np.testing.assert_array_equal(runs, [2, 3])
        np.testing.assert_array_equal(primaries, [PRIMARIES, PRIMARIES])
        self.assertEqual(means.shape, (2, bks.N_CHANNELS))
        for i, values in enumerate(RUNS.values()):
            dose = np.array([d for d, _ in values])
            error = np.array([e for _, e in values])
            n = len(values)
            np.testing.assert_allclose(means[i], dose.mean(), rtol=1e-9)
            np.testing.assert_allclose(errors[i], math.sqrt(np.sum(error**2)) / n, rtol=1e-9)
            np.testing.assert_allclose(sems[i], dose.std(ddof=1) / math.sqrt(n), rtol=1e-9)

    def test_summary_records_each_energys_runs_and_primaries(self):
        make_sweep(self.root, {"001MeVProton": 2, "002MeVProton": 3})
        energies, means, errors, sems, runs, primaries = bks.aggregate(self.root, 12)
        output = self.root / "summary.csv"
        bks.write_summary(output, energies, means, errors, sems, runs, primaries)
        with output.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(list(rows[0]), bks.FIELDNAMES)
        self.assertEqual(len(rows), 2 * bks.N_CHANNELS)
        per_energy = {(r["energy_MeV"], r["runs"], r["total_primaries"]) for r in rows}
        self.assertEqual(per_energy, {("1", "2", str(2 * PRIMARIES)),
                                      ("2", "3", str(3 * PRIMARIES))})

    def test_scalar_counts_still_accepted_by_write_summary(self):
        make_sweep(self.root, {"001MeVProton": 2, "002MeVProton": 3})
        energies, means, errors, sems, _, _ = bks.aggregate(self.root, 0)
        output = self.root / "summary.csv"
        bks.write_summary(output, energies, means, errors, sems, 5, PRIMARIES)
        with output.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual({(r["runs"], r["total_primaries"]) for r in rows},
                         {("5", str(5 * PRIMARIES))})

    def test_fixed_count_still_overrides_the_manifest(self):
        make_sweep(self.root, {"001MeVProton": 2, "002MeVProton": 3})
        with self.assertRaisesRegex(RuntimeError, "found 3 result files, expected 2"):
            bks.aggregate(self.root, 0, runs_per_energy=2, primaries_per_run=PRIMARIES)

    def test_missing_replica_is_refused(self):
        make_sweep(self.root, {"001MeVProton": 2, "002MeVProton": 4})
        with self.assertRaisesRegex(RuntimeError, "found 3 result files, expected 4"):
            bks.aggregate(self.root, 0)

    def test_wrong_entry_count_is_refused(self):
        make_sweep(self.root, {"001MeVProton": 2, "002MeVProton": 3})
        with self.assertRaisesRegex(RuntimeError, "unexpected entry count"):
            bks.aggregate(self.root, 0, primaries_per_run=2 * PRIMARIES)

    def test_manifest_without_counts_needs_fixed_values(self):
        make_sweep(self.root, {}, with_counts=False)
        with self.assertRaisesRegex(RuntimeError, "no 'array_tasks'"):
            bks.aggregate(self.root, 0)
        with self.assertRaisesRegex(RuntimeError, "no 'primaries_per_task'"):
            bks.aggregate(self.root, 0, runs_per_energy=2)
        # With both fixed the manifest columns are not needed; the file count then decides.
        with self.assertRaisesRegex(RuntimeError, "found 3 result files, expected 2"):
            bks.aggregate(self.root, 0, runs_per_energy=2, primaries_per_run=PRIMARIES)


if __name__ == "__main__":
    unittest.main()

