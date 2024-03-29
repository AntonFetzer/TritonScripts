import uproot
import awkward as ak
import numpy as np
import os
import sys

# print(f"Arguments of the script : {sys.argv[1:]=}")

# PathList = "/scratch/work/fetzera1/Gradient/Gradient/4Material/"
# Path = "/l/triton_work/Gradient/4Material/"
Path = "/scratch/work/fetzera1/MulasTest2/"
# PathList = sys.argv[1:]

# for Path in PathList:
for root, dirs, files in os.walk(os.path.abspath(Path)):
    for file in files:
        if file.endswith('.root'):
            File = os.path.join(root, file)

# Files = [f for f in os.listdir(Path + "/root/") if f.endswith('.root')]
# Files = ["5layer3-2e7proton.root"]
# Files = sys.argv[1:]
# for File in Files:

            f = uproot.open(File)
            print("Read in: " + File)

            tree = f["Detector Data 0"]

            AllKeys = tree.keys()
            keys = []

            for key in AllKeys:
                if "Edep" in key:
                    keys.append(key)

            print(keys)
            numKeys = len(keys)

            Edep = np.zeros(numKeys)
            Error = np.zeros(numKeys)

            for i in range(numKeys):
                # Edep[i] = np.sum(tree[keys[i]].array(library="np"))
                # Edep[i] = ak.sum(tree[keys[i]].array())
                Data = tree[keys[i]].array()
                Edep[i] = ak.sum(Data)
                Error[i] = ak.std(Data) * np.sqrt(len(Data))

                print(keys[i], Edep[i], Error[i])

            # np.savetxt(Path + File.split(".")[0] + ".txt", keys, Edep, Error)
            CSVFile = open(File.split(".")[0] + ".txt", 'w')
            for i in range(numKeys):
                CSVFile.writelines(keys[i] + ", " + str(Edep[i]) + ", " + str(Error[i]) + "\n")
            CSVFile.close()

# srun --mem=16G --time=05:00:00 python 100TilesMeVtoCSV.py
