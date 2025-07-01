
import sys, pkgutil, pprint, importlib
# show the search path entry you just added
print("Looking for path entry:")
pprint.pp([p for p in sys.path if 'triton_work' in p])

print("\nTrying to import GRAS …")
try:
    import GRAS
    print("GRAS found at:", GRAS.__path__)
except ImportError as e:
    print("IMPORT ERROR:", e)
    raise SystemExit

print("\nSubpackages under GRAS:")
for m in pkgutil.iter_modules(GRAS.__path__):
    print("  ", m.name)

print("\nTrying to import GRAS.Dependencies …")
try:
    deps = importlib.import_module("GRAS.Dependencies")
    print("Dependencies path:", deps.__path__)
    print("Modules inside Dependencies:")
    for m in pkgutil.iter_modules(deps.__path__):
        print("  ", m.name)
except ImportError as e:
    print("IMPORT ERROR:", e)

