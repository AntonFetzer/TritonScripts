import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def mergeHistograms(list_of_histograms):
    """
    Merges a list of GRAS histograms into a single total histogram.
    The histograms in the list must have the same collumns 'lower', 'upper', 'mean', 'value', 'error', 'entries'.
    The bins must allign.

    Args:
        list_of_histograms (list): A list of dictionaries containing the histograms to be merged.

    Returns:
        dict: A dictionary containing the merged histogram.

    Raises:
        SystemExit: If the list of histograms is empty.
        SystemExit: If the histogram does not have the expected collumns.
        SystemExit: If the histogram bins in a histogram do not match the bins in the first file.
    """

    if not list_of_histograms:
        sys.exit("ERROR !!! No histograms found")

    # NumHistograms = len(list_of_histograms)
    # print("Number of Histograms:", NumHistograms)

    # Initialize TotalHistogram for accumulation with zeros or appropriate structures
    # Since we don't have data yet, we'll initialize with None and set appropriately on first read
    TotalHistogram = {'lower': None, 'upper': None, 'mean': None, 'value': None, 'error': None, 'entries': None}

    # List to count how often a bin contained non-zero 'mean' values
    MeanCountList = None
    
    for Histogram in list_of_histograms:
        # Initialize TotalHistogram bin edges with the data of the first histogram in the list
        if TotalHistogram['lower'] is None:
            
            TotalHistogram["lower"] = Histogram["lower"]
            TotalHistogram["upper"] = Histogram["upper"]

            for key in ["mean", "value", "error", "entries"]:
                TotalHistogram[key] = np.zeros_like(Histogram[key], dtype=float)
            # Initialize MeanCountList with zero ints
            MeanCountList = np.zeros_like(Histogram['mean'], dtype=int)
        else:
            # Check if the histogram bins are aligned
            if not np.allclose(TotalHistogram['lower'], Histogram['lower']) or not np.allclose(TotalHistogram['upper'], Histogram['upper']):
                # Exit program with error message
                sys.exit(f"ERROR !!! Histogram bins in a histogram do not match the bins in the first histogram")

        # Accumulate the data
        TotalHistogram['value'] += Histogram['value'] * Histogram['entries']
        TotalHistogram['mean'] += Histogram['mean'] * Histogram['entries']
        TotalHistogram['entries'] += Histogram['entries']
        TotalHistogram['error'] += np.square(Histogram['error'] * Histogram['entries'])


    # Normalize the data and handle zero entries
    non_zero_entries = TotalHistogram['entries'] > 0
    TotalHistogram['value'][non_zero_entries] /= TotalHistogram['entries'][non_zero_entries]
    TotalHistogram['mean'][non_zero_entries] /= TotalHistogram['entries'][non_zero_entries]
    
    # Calculate the error correctly considering the weighted contributions
    TotalHistogram['error'][non_zero_entries] = np.sqrt(TotalHistogram['error'][non_zero_entries]) / TotalHistogram['entries'][non_zero_entries]

    return TotalHistogram