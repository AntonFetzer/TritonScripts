"""Combine independent measurements: ratios and inverse-variance means."""

import math


def ratioWithError(numerator, numeratorError, denominator, denominatorError):
    """Ratio of two independent measurements and its absolute error."""
    ratio = numerator / denominator
    relative = math.sqrt((numeratorError / numerator) ** 2
                         + (denominatorError / denominator) ** 2)
    return ratio, ratio * relative


def weightedMean(values, errors):
    """Inverse-variance mean, its error, and chi2 per degree of freedom about it."""
    weights = [1.0 / error ** 2 for error in errors]
    mean = sum(w * v for w, v in zip(weights, values)) / sum(weights)
    chi2 = sum(w * (v - mean) ** 2 for w, v in zip(weights, values))
    return mean, sum(weights) ** -0.5, chi2 / (len(values) - 1)
