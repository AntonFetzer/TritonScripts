import uproot
import awkward as ak
import numpy as np
import os

# "/scratch/work/fetzera1/Gradient/4MaterialGradient/Al-Pe-Al/root/",
# PathList = ["/scratch/work/fetzera1/Gradient/4MaterialGradient/Al-W-Al/root/"]

PathList = ["/scratch/work/fetzera1/Gradient/2MaterialGradient/pe-zn/root/"]

for Path in PathList:

    Files = [f for f in os.listdir(Path) if f.endswith('.root')]
    # Files = ["gradient-al-test-directions-1e7proton.root"]

    for File in Files:
        f = uproot.open(Path + File)
        print("Read in: " + Path + File)

        tree = f["Detector Data 0"]

        keys = tree.keys()

        # print(keys)

        Edep = np.zeros(len(keys))

        for i in range(len(keys)):
            # Edep[i] = np.sum(tree[keys[i]].array(library="np"))
            Edep[i] = ak.sum(tree[keys[i]].array())
            print(keys[i], Edep[i])

        np.savetxt(Path + File.split(".")[0] + ".txt", Edep)

# srun --mem=16G --time=05:00:00 python 100TilesMeVtoCSV.py
