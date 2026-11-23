"""Dependencies.MergeHistograms: the weighting of identical runs.

Values and errors are weighted by the primary events of each job. Weighting by
the entries of a bin instead drops the jobs that saw nothing in it, which once
overestimated the deep-shielding SEE rates by up to a factor of 5.
"""
import unittest

import numpy as np

from Dependencies.MergeHistograms import mergeHistograms


def job(events, mean, value, error, entries):
    return {
        'lower': np.array([1.0, 2.0]), 'upper': np.array([2.0, 3.0]), 'mean': np.array(mean),
        'value': np.array(value), 'error': np.array(error), 'entries': np.array(entries), 'events': events,
    }


class TestMergeHistograms(unittest.TestCase):

    def test_weighted_by_events_including_an_empty_job(self):
        merged = mergeHistograms([
            job(1000.0, [1.2, 2.5], [8.0, 4.0], [2.0, 4.0], [10, 5]),
            job(2000.0, [1.6, 0.0], [2.0, 0.0], [0.5, 0.0], [30, 0]),
            job(5000.0, [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0, 0]),
        ])
        # 8000 events: (1000 * 8 + 2000 * 2) / 8000 and 1000 * 4 / 8000
        np.testing.assert_allclose(merged['value'], [1.5, 0.5])
        # sqrt((1000 * 2)^2 + (2000 * 0.5)^2) / 8000 and 1000 * 4 / 8000
        np.testing.assert_allclose(merged['error'], [np.sqrt(5.0) / 8, 0.5])
        # The mean inside a bin is weighted by its entries: (10 * 1.2 + 30 * 1.6) / 40
        np.testing.assert_allclose(merged['mean'], [1.5, 2.5])
        np.testing.assert_array_equal(merged['entries'], [40, 5])
        self.assertEqual(merged['events'], 8000.0)


if __name__ == "__main__":
    unittest.main()
