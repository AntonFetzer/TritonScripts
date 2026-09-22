"""Validate installed GaAs NIEL lookups without changing the GRAS installation.

Run from the Python repository with the normal Geant4 environment loaded.
Uses Read.ReadNID.readNID for result files. Each transport check has 2000
histories and a 10-second timeout, is single-threaded, and is deliberately
a focused detector diagnostic, NOT a production HUS beam simulation.
Temporary macros, results and logs are removed after saving evidence.
After the permanent 05.02.01 launcher is installed, use mode "installed" to
validate normal gras. Historical modes use the preserved original executable.
NIEL_VALIDATION_OUTPUT selects a separate evidence filename within ROOT.
"""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from Read.ReadNID import readNID

ROOT = Path("/scratch/work/fetzera1/GRAS/CRDS/GRAS-NIEL-Bug-Analysis")
PROD = Path("/scratch/work/fetzera1/GRAS/CRDS/CRDS-HUS/PDD-Electron-CRDS-20mmField")
MAIN = Path("/scratch/work/fetzera1/GRAS/CRDS/CRDS-main.mac")


def validate():
    aliases = (PROD / "A.mac").read_text().split("/control/execute ../../CRDS-main.mac")[0]
    main = "\n".join(line for line in MAIN.read_text().splitlines()
                     if not line.startswith("/control/doifBatch"))
    extra = ""
    for suffix, volume in [("Active", "LED_active_0_PV"), ("Substrate", "LED_die_0_PV")]:
        name = "summers" + suffix
        extra += (
            f"\n/gras/analysis/nid/addModule {name}"
            f"\n/gras/analysis/nid/{name}/addVolume {volume}"
            f"\n/gras/analysis/nid/{name}/setCurveSet GaAs"
            f"\n/gras/analysis/nid/{name}/bookHistos false"
            f"\n/gras/analysis/nid/{name}/bookTuples false\n"
        )
    destination = ROOT / os.environ.get("NIEL_VALIDATION_OUTPUT", "transport-validation.json")
    evidence = json.loads(destination.read_text()) if destination.exists() else {}
    for energy in [5.0, 90.0]:
        for mode in os.environ.get("NIEL_VALIDATION_MODES", "native,observe,fixed,fix-only").split(","):
            tag = f"{energy:g}MeV-{mode}"
            with tempfile.TemporaryDirectory(prefix=".validation-", dir=ROOT) as directory:
                folder = Path(directory)
                (folder / "Res").mkdir()
                # The near-die source makes a short, useful scoring check possible.
                source = (
                    "\n/gps/source/clear\n/gps/source/add 1\n/gps/particle e-"
                    "\n/gps/number 1\n/gps/pos/type Plane\n/gps/pos/shape Rectangle"
                    "\n/gps/pos/centre 5 0 1.00 mm"
                    "\n/gps/pos/halfx 0.15 mm\n/gps/pos/halfy 0.15 mm"
                    "\n/gps/ang/type planar\n/gps/direction 0 0 -1"
                    f"\n/gps/ene/type Mono\n/gps/ene/mono {energy} MeV"
                    f"\n/gras/histo/fileName {folder}/Res/TID"
                    "\n/gras/event/printModulo 2000\n/run/beamOn 2000\n"
                )
                macro = folder / "A.mac"
                switch = ""
                if mode == "macro-only":
                    switch = ("\n/gras/analysis/nid/nidActive/setCurveSet GaAs"
                              "\n/gras/analysis/nid/nidSubstrate/setCurveSet GaAs\n")
                macro.write_text(aliases + main + extra + switch + source)
                environment = dict(os.environ, LD_PRELOAD=str(ROOT / "libniel_audit.so"),
                                   GRAS_NIEL_FIX="1" if mode == "fixed" else "0",
                                   OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
                if mode in ("native", "macro-only", "installed"):
                    environment.pop("LD_PRELOAD", None)
                elif mode == "fix-only":
                    environment["LD_PRELOAD"] = str(ROOT / "libniel_fix.so")
                original = Path("/home/fetzera1/software/gras/05.02.01/install/bin/gras-05.02.01-original")
                executable = str(original) if original.exists() and mode != "installed" else "gras"
                try:
                    run = subprocess.run([executable, str(macro)], cwd=PROD, env=environment,
                                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                         text=True, timeout=10)
                except subprocess.TimeoutExpired as exc:
                    output = exc.stdout or b""
                    if isinstance(output, bytes):
                        output = output.decode(errors="replace")
                    evidence[tag] = {"timed_out_seconds": 10,
                                     "last_lines": output.splitlines()[-40:],
                                     "audit": [line for line in output.splitlines()
                                               if line.startswith(("NIEL_AUDIT", "NIEL_FIX"))]}
                    print(tag, "TIMED OUT", flush=True)
                    continue
                lines = run.stdout.splitlines()
                csv_files = list((folder / "Res").glob("*.csv"))
                if run.returncode or len(csv_files) != 1:
                    raise RuntimeError(f"{tag}: exit {run.returncode}; files {csv_files}; "
                                       + "\n".join(lines[-30:]))
                modules = readNID(csv_files[0])
                expected = {"nidActive", "nidSubstrate", "summersActive", "summersSubstrate"}
                if set(modules) != expected or any(v["entries"] != 2000 for v in modules.values()):
                    raise RuntimeError(f"Incomplete validation: {tag}: {modules}")
                evidence[tag] = {
                    "energy_MeV": energy,
                    "description": "0.30 mm square plane at (5,0,1.00) mm; diagnostic only",
                    "exit_code": run.returncode,
                    "executable": executable,
                    "modules": modules,
                    "audit": [line for line in lines if line.startswith(("NIEL_AUDIT", "NIEL_FIX"))],
                    "normalisation": [line for line in lines if "Final normalisation factor:" in line],
                    "warnings": [line for line in lines if
                                 ("ERROR" in line or "Exception" in line or
                                  "COMMAND NOT FOUND" in line) and
                                 not line.startswith("#")],
                }
                print(tag, json.dumps(modules), flush=True)
    destination = ROOT / os.environ.get("NIEL_VALIDATION_OUTPUT", "transport-validation.json")
    destination.write_text(json.dumps(evidence, indent=2) + "\n")
    print(f"Wrote {destination}")


if __name__ == "__main__":
    validate()
