import uproot
import awkward as ak
# import numpy as np
import timeit
import os
from concurrent.futures import ThreadPoolExecutor
import multiprocessing

print("Root Version:", uproot.__version__)

executor = uproot.ThreadPoolExecutor(num_workers=os.cpu_count())  # threads

print("Number of Threads", os.cpu_count())

Path = "/l/RootScripts/root/"
# Path = "/scratch/work/fetzera1/Gradient/4MaterialGradient/Al-Pe-Al/root/"

# Files = [f for f in os.listdir(Path) if f.endswith('.root')]
Files = ["gradientalal22e9electron.root"]

start = timeit.default_timer()

for File in Files:
    f = uproot.open(Path + File)
    print("Read in: " + Path + File)

    tree = f["Detector Data 0"]

    keys = tree.keys()

    # print(keys)

    # Edep = np.zeros(len(keys))

    Edep = []

    # Edep = Parallel(n_jobs=40)(delayed(np.sum(tree[keys[i]].array(library="np"))) for i in range(len(keys)))

    # Data = tree.arrays(executor=executor)

    # print(type(Data))
    for i in range(len(keys)):
        Edep.append(ak.sum(tree[keys[i]].array()))
        print(keys[i], Edep[i])

    # np.savetxt(Path + File.split(".")[0] + ".txt", Edep)

stop = timeit.default_timer()

print('Time: ', stop - start)

# srun --mem=16G --time=00:30:00 python ParalelTest.py
