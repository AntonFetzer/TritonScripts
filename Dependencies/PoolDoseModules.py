"""Pool repeated GRAS single-volume dose modules without changing their units."""
import numpy as np
from Read.ReadDose import readDoseModules

def poolDoseModules(files, module_prefix="doseLayer-"):
    """Entry-weight module doses; return pooled values and replica Birge diagnostics.

    Preserves source units. Requires matching modules/units, positive entries,
    and finite nonnegative dose/errors. Birge > 2 is flagged, never concealed.
    The caller must investigate flags before treating the dataset as validated.
    """
    data = [readDoseModules(p, module_prefix) for p in files]
    if not data or not data[0]:
        raise ValueError("No dose modules")
    names = set(data[0])
    if any(set(d) != names for d in data):
        raise ValueError("Module mismatch")
    result = {}
    for name in sorted(names):
        items = [d[name] for d in data]
        units = {d["unit"] for d in items}
        if len(units) != 1:
            raise ValueError("Unit mismatch")
        n = np.array([d["entries"] for d in items], dtype=float)
        y = np.array([d["dose"] for d in items])
        s = np.array([d["error"] for d in items])
        nz = np.array([d["non-zeros"] for d in items])
        if not np.all(np.isfinite([n,y,s,nz])) or np.any(n<=0) or np.any(y<0) or np.any(s<0) or np.any(nz<0) or np.any(nz>n):
            raise ValueError("Invalid dose data")
        mean = np.sum(n*y)/n.sum()
        internal = np.sum((n*s)**2)
        external = np.sum((n*(y-mean))**2)
        birge = np.sqrt(external/internal) if internal else (0. if external == 0 else float("inf"))
        result[name] = dict(dose=mean, error=np.sqrt(internal)/n.sum(),
                            entries=int(n.sum()), unit=units.pop(),
                            birge_ratio=birge, inconsistent=bool(birge>2),
                            **{"non-zeros":int(nz.sum())})
    return result
