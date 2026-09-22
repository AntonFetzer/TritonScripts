"""Delivered TID and DDD from per-fluence coefficients and local fluence maps.

A GRAS coefficient is per incident particle/cm2. Where a campaign delivered
several exposures with fields that were not uniform over the instrument, the
delivered quantity at each part is the sum over exposures of the coefficient
times the fluence AT THAT PART in that exposure:

    delivered = sum_e  C(run_e) * phi_e(part)

The local fluence phi_e(part) is sampled from the accelerator-log isocentre
map of exposure e with UppsalaAcceleratorLogs.radex_isocentre_outputs
.bilinear_sample, the sampler behind the RadEx per-RadFET fluence table. The
field-average fluence, delivered protons over requested area, is carried
alongside for comparison.

Coefficients from one run are fully correlated between the exposures that use
that run, so the statistical error is combined per run first:

    error = sqrt( sum_r ( sigma_C(r) * sum_{e using r} phi_e )^2 )

Board alignment on the field was not recorded. The sensitivity to it is the
spread of the delivered value when the whole board is shifted by +-d in x and
in y, all exposures together.
"""

import math
from pathlib import Path

from Read.ReadCRDSFluence import fieldAverageFluence, sampleFluence


def exposureFluences(delivery, shift=(0.0, 0.0)):
    """Local and field-average fluence per exposure and part.

    Args:
        delivery (dict): a campaign's ``delivery`` block.
        shift (tuple): board displacement (dx, dy) in mm applied to every part.

    Returns:
        list[dict]: per exposure, ``{"label", "fieldAverage", "local": {part: phi}}``.
    """
    maps = Path(delivery["maps"])
    exposures = []
    for exposure in delivery["exposures"]:
        local = {part: sampleFluence(maps / exposure["map"] / "heatmap_isocentre.npz",
                                     x + shift[0], y + shift[1])
                 for part, (x, y) in delivery["parts"].items()}
        exposures.append({
            "label": exposure["label"],
            "fieldAverage": fieldAverageFluence(maps / exposure["map"] / "summary.json"),
            "local": local,
        })
    return exposures


def combine(coefficients, runOf, fluences):
    """Delivered value and statistical error for one tally.

    Args:
        coefficients (dict): run label -> (coefficient, error) for this tally.
        runOf (list[str]): per exposure, the run whose coefficient it uses.
        fluences (list[float]): per exposure, the fluence at this tally's part.

    Returns:
        tuple: (delivered, error), errors combined per run as described above.
    """
    if len(runOf) != len(fluences):
        raise ValueError("one run and one fluence per exposure are required")
    perRun = {}
    for run, fluence in zip(runOf, fluences):
        perRun[run] = perRun.get(run, 0.0) + fluence
    delivered = sum(coefficients[run][0] * total for run, total in perRun.items())
    error = math.sqrt(sum((coefficients[run][1] * total) ** 2
                          for run, total in perRun.items()))
    if not (math.isfinite(delivered) and math.isfinite(error)):
        raise ValueError("non-finite delivered value")
    return delivered, error


def deliveredRecords(delivery, coefficients):
    """Delivered TID and DDD per tally, with the comparisons that qualify it.

    Args:
        delivery (dict): a campaign's ``delivery`` block.
        coefficients (dict): (quantity, tally) -> {run label: (value, error)}.

    Returns:
        list[dict]: per tally, the delivered value from the local fluence, its
            statistical error, the per-exposure contributions, the value from
            the field-average fluence, and the range over board shifts of
            +-alignmentShift mm.
    """
    nominal = exposureFluences(delivery)
    step = delivery["alignmentShift_mm"]
    shifted = [exposureFluences(delivery, shift)
               for shift in ((step, 0.0), (-step, 0.0), (0.0, step), (0.0, -step))]

    records = []
    for (quantity, tally), runs in coefficients.items():
        part = delivery["partOfQuantity"][quantity]
        runOf = [exposure["coefficientRun"][part] for exposure in delivery["exposures"]]
        local = [exposure["local"][part] for exposure in nominal]
        delivered, error = combine(runs, runOf, local)
        average, _ = combine(runs, runOf, [exposure["fieldAverage"] for exposure in nominal])
        aligned = [combine(runs, runOf, [exposure["local"][part] for exposure in fluences])[0]
                   for fluences in shifted]
        records.append({
            "quantity": quantity,
            "tally": tally,
            "part": part,
            "exposures": [{
                "label": exposure["label"],
                "run": run,
                "coefficient": runs[run][0],
                "coefficientError": runs[run][1],
                "localFluence": exposure["local"][part],
                "fieldAverageFluence": exposure["fieldAverage"],
                "delivered": runs[run][0] * exposure["local"][part],
            } for exposure, run in zip(nominal, runOf)],
            "delivered": delivered,
            "error": error,
            "fieldAverageDelivered": average,
            "alignmentLow": min(aligned),
            "alignmentHigh": max(aligned),
        })
    return records
