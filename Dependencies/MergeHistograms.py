import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def mergeHistograms(list_of_histograms):
    """
    Merges a list of GRAS histograms of identical runs using proper error
    propagation.

    Each job's normalized values are weighted by the job's number of primary
    events ('events' key, parsed from the GRAS csv NumOfEvt field), falling
    back to the job's TOTAL entry count for histograms without it.  Weighting
    by per-bin entries instead would drop jobs that saw no event in a bin,
    which inflates rare-event bins by up to N_jobs/N_jobs_with_hits (this
    overestimated the deep-shielding SEE rates by up to x5).  A completed job
    with an all-zero histogram is a valid zero measurement and correctly
    pulls the average down.

    The per-bin 'mean' stays entry-weighted: it locates the average LET of
    the recorded events inside the bin, which is only defined given a hit.

    Args:
        list_of_histograms (list): List of histograms with aligned bins

    Returns:
        dict: Merged histogram with proper error propagation
    """
    if not list_of_histograms:
        sys.exit("ERROR !!! No histograms found")

    # Validate that all histograms have the required keys
    required_keys = ['lower', 'upper', 'mean', 'value', 'error', 'entries']
    for hist in list_of_histograms:
        if not all(key in hist for key in required_keys):
            sys.exit("ERROR !!! Missing required histogram keys")

    # Check that all histograms have the same length (before allclose, which
    # raises on mismatched shapes instead of reporting cleanly)
    bin_count = len(list_of_histograms[0]['lower'])
    for hist in list_of_histograms[1:]:
        if len(hist['lower']) != bin_count:
            sys.exit("ERROR !!! Mismatched histogram lengths")

    # Check that all histograms have the same binning
    for hist in list_of_histograms[1:]:
        if not np.allclose(list_of_histograms[0]['lower'], hist['lower']) or \
           not np.allclose(list_of_histograms[0]['upper'], hist['upper']):
            sys.exit("ERROR !!! Mismatched histogram bins")

    # Check that all histograms have the correct data types
    # Floats for 'bins', 'value', 'error', 'mean' and integers for 'entries'
    for hist in list_of_histograms:
        for key in ['lower', 'upper', 'mean', 'value', 'error']:
            if not np.issubdtype(np.array(hist[key]).dtype, np.floating):
                sys.exit(f"ERROR !!! Histogram key '{key}' must be of float type")
        if not np.issubdtype(np.array(hist['entries']).dtype, np.integer):
            sys.exit("ERROR !!! Histogram key 'entries' must be of integer type")


    TotalHistogram = {
        'lower': list_of_histograms[0]['lower'].copy(),
        'upper': list_of_histograms[0]['upper'].copy(),
        'value': np.zeros_like(list_of_histograms[0]['value']),
        'mean': np.zeros_like(list_of_histograms[0]['mean']), 
        'entries': np.zeros_like(list_of_histograms[0]['entries']),
        'error': np.zeros_like(list_of_histograms[0]['error'])
    }

    # First pass: accumulate weighted sums.  Values and errors are weighted
    # by the job's total entries; the in-bin mean by the bin's entries.
    total_weight = 0.0
    weights = []
    for hist in list_of_histograms:
        if not np.allclose(TotalHistogram['lower'], hist['lower']) or \
           not np.allclose(TotalHistogram['upper'], hist['upper']):
            sys.exit("ERROR !!! Mismatched histogram bins")

        entries = hist['entries']
        weight = float(hist.get('events', 0)) or float(np.sum(entries))
        weights.append(weight)
        total_weight += weight
        TotalHistogram['entries'] += entries
        TotalHistogram['value'] = TotalHistogram['value'] + hist['value'] * weight
        TotalHistogram['mean'] += hist['mean'] * entries
        TotalHistogram['error'] = TotalHistogram['error'] + (hist['error'] * weight) ** 2

    # Second pass: normalize.  total_weight can only be zero if every job is
    # empty AND carries no event count; the all-zero histogram is returned.
    if total_weight > 0:
        TotalHistogram['value'] /= total_weight
        InternalErrSq = TotalHistogram['error'].copy()  # sum of (w_i * sigma_i)^2
        TotalHistogram['error'] = np.sqrt(TotalHistogram['error']) / total_weight
    else:
        print("WARNING: merging histograms with zero total weight")
    mask = TotalHistogram['entries'] > 0
    TotalHistogram['mean'][mask] /= TotalHistogram['entries'][mask]

    # Consistency diagnostic (Birge ratio): for jobs that differ only by their
    # random seed, the observed scatter of the per-job values around the merged
    # value is fully explained by the quoted statistical errors, so
    # external/internal ~ 1 per bin.  A ratio well above 1 in a well-populated
    # bin means the jobs genuinely differ (cuts, source sampling, build,
    # biasing, or a corrupt csv) -- that is a bug signal, not extra statistics,
    # so halt and let the user look instead of silently merging.
    # Only bins with enough entries are checked; in sparser bins the ratio
    # itself is too noisy to judge (spread ~ sqrt(2/hits)).
    MinEntriesForBirgeCheck = 100
    BirgeThreshold = 2.0
    if total_weight > 0 and len(list_of_histograms) > 1:
        ExternalErrSq = np.zeros_like(TotalHistogram['value'])
        for hist, weight in zip(list_of_histograms, weights):
            ExternalErrSq += weight ** 2 * (hist['value'] - TotalHistogram['value']) ** 2
        Birge = np.sqrt(np.divide(ExternalErrSq, InternalErrSq,
                                  out=np.zeros_like(InternalErrSq), where=InternalErrSq > 0))
        Bad = (TotalHistogram['entries'] >= MinEntriesForBirgeCheck) & (Birge > BirgeThreshold)
        if np.any(Bad):
            print("ERROR !!! Inconsistent jobs: the between-job scatter exceeds the quoted statistical errors")
            for b in np.flatnonzero(Bad):
                print(f"  bin {b} [{TotalHistogram['lower'][b]:.4g}, {TotalHistogram['upper'][b]:.4g}]:"
                      f" Birge ratio {Birge[b]:.2f} with {TotalHistogram['entries'][b]} entries")
            print("  The jobs differ by more than their random seed, or a csv file is corrupt.")
            try:
                input("  Press Enter to continue with the merge anyway...")
            except EOFError:
                print("  (non-interactive run, continuing)")

    # Total number of primary events, where the inputs provide it
    TotalHistogram['events'] = float(sum(hist.get('events', 0) for hist in list_of_histograms))

    return TotalHistogram