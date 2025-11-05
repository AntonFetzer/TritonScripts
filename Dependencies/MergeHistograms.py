import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def mergeHistograms(list_of_histograms):
    """
    Merges a list of GRAS histograms using proper error propagation.
    
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

    # Check that all histograms have the same binning
    for hist in list_of_histograms[1:]:
        if not np.allclose(list_of_histograms[0]['lower'], hist['lower']) or \
           not np.allclose(list_of_histograms[0]['upper'], hist['upper']):
            sys.exit("ERROR !!! Mismatched histogram bins")

    # Check that all histograms have the same length
    bin_count = len(list_of_histograms[0]['lower'])
    for hist in list_of_histograms[1:]:
        if len(hist['lower']) != bin_count:
            sys.exit("ERROR !!! Mismatched histogram lengths")

    # Check that all histograms have the correct data types
    # Floats for 'bins', 'value', 'error', 'mean' and integers for 'entries'
    for hist in list_of_histograms:
        for key in ['lower', 'upper', 'mean', 'value', 'error']:
            if not np.issubdtype(np.array(hist[key]).dtype, np.floating):
                sys.exit(f"ERROR !!! Histogram key '{key}' must be of float type")
        if not np.issubdtype(np.array(hist['entries']).dtype, np.integer):
            sys.exit("ERROR !!! Histogram key 'entries' must be of integer type")

    # Check that all histograms have non zero total value
    for hist in list_of_histograms:
        if np.sum(hist['value']) == 0:
            sys.exit("ERROR !!! One of the histograms has zero total value")

    TotalHistogram = {
        'lower': list_of_histograms[0]['lower'].copy(),
        'upper': list_of_histograms[0]['upper'].copy(),
        'value': np.zeros_like(list_of_histograms[0]['value']),
        'mean': np.zeros_like(list_of_histograms[0]['mean']), 
        'entries': np.zeros_like(list_of_histograms[0]['entries']),
        'error': np.zeros_like(list_of_histograms[0]['error'])
    }

    # First pass: accumulate weighted sums
    for hist in list_of_histograms:
        if not np.allclose(TotalHistogram['lower'], hist['lower']) or \
           not np.allclose(TotalHistogram['upper'], hist['upper']):
            sys.exit("ERROR !!! Mismatched histogram bins")
            
        entries = hist['entries']
        TotalHistogram['entries'] += entries
        TotalHistogram['value'] += hist['value'] * entries  
        TotalHistogram['mean'] += hist['mean'] * entries
        TotalHistogram['error'] += (hist['error'] * entries) ** 2

    # Second pass: normalize
    mask = TotalHistogram['entries'] > 0
    TotalHistogram['value'][mask] /= TotalHistogram['entries'][mask]
    TotalHistogram['mean'][mask] /= TotalHistogram['entries'][mask]
    TotalHistogram['error'][mask] = np.sqrt(TotalHistogram['error'][mask]) / TotalHistogram['entries'][mask]

    return TotalHistogram