"""Read the settings of one GRAS run from its own macros and logs.

A run directory holds A.mac, which sets aliases and executes the study's main
macro, and GRAS runs inside that directory, so every /control/execute path
resolves against it. runSettings follows that chain, substituting {aliases} as
Geant4 does, and keeps the last argument given to every command. The analysis
tools read the source spectrum, the source area and the normalisation from it
instead of hard-coding them per run, so pointing a tool at another run
directory is enough.

Loops, conditionals and interactive branches (/control/loop, /control/foreach,
/control/if*, /control/doif*) are not followed: in the RadEx main macros they
only add volumes to regions, visualise or start the run, and none of them sets
the source or the normalisation.
"""

from pathlib import Path
import re

#: Unit factors to MeV and to mm.
ENERGY_UNITS = {"eV": 1e-6, "keV": 1e-3, "MeV": 1.0, "GeV": 1e3}
LENGTH_UNITS = {"um": 1e-3, "mm": 1.0, "cm": 10.0, "m": 1e3}
#: Slurm log names: run.sh-job<array job>-<task>.log (current convention) and
#: <name>-<array job>_<task>.out/.err (older scripts).
LOG_PATTERNS = (re.compile(r"job(\d+)-\d+\.log$"), re.compile(r"-(\d+)_\d+\.(?:out|err)$"))

_ALIAS = re.compile(r"\{(\w+)\}")


def _substitute(text, aliases):
    return _ALIAS.sub(lambda m: aliases.get(m.group(1), m.group(0)), text)


def runSettings(runPath):
    """Aliases and commands in effect after a run's macro chain.

    Args:
        runPath (str or Path): run directory holding A.mac.

    Returns:
        dict: 'aliases' (name -> value), 'commands' (command -> its last
            argument string, aliases substituted) and 'macros' (every macro
            read, in order).

    Raises:
        FileNotFoundError: A.mac or an executed macro does not exist.
    """
    run = Path(runPath)
    aliases, commands, macros = {}, {}, []

    def walk(macro):
        if not macro.is_file():
            raise FileNotFoundError(f"{macro} (executed from {macros[-1] if macros else run})")
        macros.append(macro)
        for raw in macro.read_text(encoding="utf-8").splitlines():
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            command, _, argument = _substitute(line, aliases).partition(" ")
            argument = argument.strip()
            if command == "/control/alias":
                name, _, value = argument.partition(" ")
                aliases[name] = value.strip().strip('"')
            elif command == "/control/multiply":
                name, a, b = argument.split()
                aliases[name] = repr(float(a) * float(b))
            elif command == "/control/execute":
                walk((run / argument).resolve())
            else:
                commands[command] = argument

    walk(run / "A.mac")
    return {"aliases": aliases, "commands": commands, "macros": macros}


def _quantity(argument, units, what):
    value, _, unit = argument.partition(" ")
    unit = unit.strip()
    if unit not in units:
        raise ValueError(f"unknown {what} unit {unit!r} in {argument!r}")
    return float(value) * units[unit]


def spectrum(runPath, settings=None):
    """Source energy distribution of a single-source GPS run.

    Returns:
        dict: 'type' (GPS ene/type, 'Mono' when unset), 'centre_MeV' and
            'sigma_MeV' (0 when unset).
    """
    commands = (settings or runSettings(runPath))["commands"]
    if "/gps/source/add" in commands:
        raise ValueError(f"{runPath}: multi-source GPS; read it with Read.ReadGPSSources")
    if "/gps/ene/mono" not in commands:
        raise ValueError(f"{runPath}: no /gps/ene/mono in its macros")
    sigma = commands.get("/gps/ene/sigma")
    return {
        "type": commands.get("/gps/ene/type", "Mono"),
        "centre_MeV": _quantity(commands["/gps/ene/mono"], ENERGY_UNITS, "energy"),
        "sigma_MeV": _quantity(sigma, ENERGY_UNITS, "energy") if sigma else 0.0,
    }


def sourceAreaCm2(runPath, settings=None):
    """Area of a planar rectangular or square GPS source [cm2]."""
    commands = (settings or runSettings(runPath))["commands"]
    shape = commands.get("/gps/pos/shape")
    if commands.get("/gps/pos/type") != "Plane" or shape not in ("Rectangle", "Square"):
        raise ValueError(f"{runPath}: source is not a planar rectangle or square")
    halfx = _quantity(commands["/gps/pos/halfx"], LENGTH_UNITS, "length")
    halfy = halfx if shape == "Square" else _quantity(commands["/gps/pos/halfy"], LENGTH_UNITS, "length")
    return 4.0 * halfx * halfy / 100.0


def normalisationFluence(runPath, settings=None):
    """Incident fluence GRAS normalised the run to [1/cm2], from setSourceFluence."""
    commands = (settings or runSettings(runPath))["commands"]
    argument = commands.get("/gras/analysis/setSourceFluence")
    if argument is None:
        raise ValueError(f"{runPath}: no /gras/analysis/setSourceFluence in its macros")
    value, _, unit = argument.partition(" ")
    if unit.strip() != "1/cm2":
        raise ValueError(f"{runPath}: source fluence unit {unit.strip()!r}, expected 1/cm2")
    return float(value)


def slurmJobs(runPath):
    """Sorted Slurm array-job IDs named by the run's log/ files; empty if it has none."""
    jobs = set()
    for path in (Path(runPath) / "log").glob("*"):
        for pattern in LOG_PATTERNS:
            match = pattern.search(path.name)
            if match:
                jobs.add(int(match.group(1)))
    return sorted(jobs)
