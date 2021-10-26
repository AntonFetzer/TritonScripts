import uproot
import awkward as ak
import numpy as np
import os


PathList = ["/scratch/work/fetzera1/Permutations/3Layer1-5gcm2/"]
# PathList = ["/home/anton/Desktop/triton_work/TEST/100Materials/root/"]

for Path in PathList:

    Files = [f for f in os.listdir(Path) if f.endswith('.root')]
    # Files = ["100materials2e9electron.root"]

    N = 100

    for File in Files:
        f = uproot.open(Path + File)
        print("Read in: " + Path + File)

        tree = f["Detector Data 0"]

        keys = tree.keys()

        # print(keys)
        numKeys = len(keys)

        Edep = np.zeros(numKeys)
        Error = np.zeros(numKeys)

        for i in range(numKeys):
            # Edep[i] = np.sum(tree[keys[i]].array(library="np"))
            # Edep[i] = ak.sum(tree[keys[i]].array())
            Data = tree[keys[i]].array()
            Edep[i] = ak.sum(Data)

            SampleLen = int(len(Data) / N)
            Sums = np.zeros(N)

            for j in range(N):
                Sums[j] = sum(Data[SampleLen * j:SampleLen * (j + 1) - 1])
                # print("From:", SampleLen*j)
                # print("To:", SampleLen*(j+1)-1)
                # print("Sum:", Samples[j])

            Error[i] = np.std(Sums) * np.sqrt(N)

            print(keys[i], Edep[i], Error[i])

        #np.savetxt(Path + File.split(".")[0] + ".txt", keys, Edep, Error)
        CSVFile = open(Path + File.split(".")[0] + ".txt", 'w')
        for i in range(numKeys):
            CSVFile.writelines(keys[i] + ", " + str(Edep[i]) + ", " + str(Error[i]) + "\n")
        CSVFile.close()

# srun --mem=16G --time=05:00:00 python 100TilesMeVtoCSV.py
