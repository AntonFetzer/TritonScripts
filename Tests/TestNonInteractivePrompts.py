"""Regression test for the "Press Enter" prompts of the shared readers.

mergeHistograms, totalDose and totalFluenceHistos ask for Enter when the
merged files are inconsistent or missing. Run over ssh or in a pipeline, stdin
can be an open pipe that never delivers EOF, and input() then blocks forever
(a fluence plot hung this way). The prompts must therefore only appear when
stdin is a terminal.

Run from the repository root:

    python3 -m unittest discover -s Tests -p "Test*.py"
"""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.MergeHistograms import mergeHistograms
from Dependencies.TotalFluenceHistograms import totalFluenceHistos


def histogram(value):
    """One-bin histogram with a small quoted error and many entries."""
    return {
        'lower': np.array([1.0]), 'upper': np.array([2.0]), 'mean': np.array([1.5]),
        'value': np.array([value]), 'error': np.array([0.01]),
        'entries': np.array([1000]), 'events': 1000.0,
    }


class TestNonInteractivePrompts(unittest.TestCase):

    def mergeInconsistent(self, isatty):
        """Merge two jobs whose scatter far exceeds their errors (Birge ratio >> 2)."""
        with mock.patch("sys.stdin") as stdin, \
             mock.patch("builtins.input") as prompt, \
             mock.patch("builtins.print"):
            stdin.isatty.return_value = isatty
            merged = mergeHistograms([histogram(1.0), histogram(2.0)])
        return merged, prompt

    def test_merge_does_not_prompt_without_terminal(self):
        merged, prompt = self.mergeInconsistent(isatty=False)
        prompt.assert_not_called()
        self.assertAlmostEqual(merged['value'][0], 1.5)

    def test_merge_prompts_on_terminal(self):
        _, prompt = self.mergeInconsistent(isatty=True)
        prompt.assert_called_once()

    def test_missing_fluence_files_do_not_prompt_without_terminal(self):
        with tempfile.TemporaryDirectory() as empty, \
             mock.patch("sys.stdin") as stdin, \
             mock.patch("builtins.input") as prompt, \
             mock.patch("builtins.print"):
            stdin.isatty.return_value = False
            self.assertIsNone(totalFluenceHistos(empty))
        prompt.assert_not_called()


if __name__ == "__main__":
    unittest.main()
