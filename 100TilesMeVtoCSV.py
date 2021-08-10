import uproot
import numpy as np
import os

Path = "/home/anton/Desktop/triton_work/Gradient/4MaterialGradient/"
# Path = "/scratch/work/fetzera1/Gradient/2MaterialGradient/ta-al/root/"

Files = [f for f in os.listdir(Path) if f.endswith('.root')]
#Files = ["gradient-al-pe-al-pe-1e8electron.root"]

for File in Files:
    f = uproot.open(Path + File)
    print("Read in: " + Path + File)

    tree = f["Detector Data 0"]

    keys = tree.keys()

    print(keys)

    Edep = np.zeros(len(keys))

    for i in range(len(keys)):
        Edep[i] = np.sum(tree[keys[i]].array(library="np"))

    np.savetxt(Path + File.split(".")[0] + ".txt", Edep)


# srun --mem=10G --time=00:15:00 python 100TilesMeVtoCSV.py
