from ReadG4root import readG4root
import numpy as np
import glob
import matplotlib.pyplot as plt

# ---------------------- Dont Modify -------- Dont break compatability -----------------------

# ---------------------- Dont Modify -------- Dont break compatability -----------------------

# ---------------------- Dont Modify -------- Dont break compatability -----------------------

Folder = "/scratch/work/fetzera1/MULASS/MulasCosmicProton1e8Full/"

# Find all files ending with .root
FileList = glob.glob(Folder + "root/*.root")

for file in FileList:
    Data = readG4root(file)

    Length = len(Data[0])
    print("Number of Data points", "{:,}".format(Length))
    Len = int(Length / 100)

    Samples = np.zeros((3, 100))

    for i in range(100):
        Samples[0, i] = sum(Data[0, Len * i:Len * (i + 1) - 1])
        Samples[1, i] = sum(Data[1, Len * i:Len * (i + 1) - 1])
        Samples[2, i] = sum(Data[2, Len * i:Len * (i + 1) - 1])
        # print("From:", Len*i)
        # print("To:", Len*(i+1)-1)

    Results = [['' for i in range(5)] for j in range(4)]

    print("Results", Results)

    # 0: Labels
    Results[0] = ["Branch", "Total MeV", "Std Dev MeV", "Std Dev %", "Num of entries per chunk"]
    Results[1][0] = "GunEnergy"
    Results[2][0] = "Edep"
    Results[3][0] = "Esec"

    for i in range(3):
        # 1: Total MeV
        Results[i + 1][1] = sum(Data[i])
        # 2: Standard Deviation MeV
        Results[i + 1][2] = np.std(Samples[i, :])
        # 3: Standard Deviation %
        Results[i + 1][3] = 100 * Results[i + 1][2] / Results[i + 1][1]
        # 4:Num of entries per chunk
        Results[i + 1][4] = Length

    print("Results[0]", Results[0])
    print("Results[1]", Results[1])
    print("Results[2]", Results[2])
    print("Results[3]", Results[3])

    CSVfile = open(file.split(".")[0] + ".csv", 'w')
    for i in range(4):
        CSVfile.write(','.join(str(e) for e in Results[i]) + "\n")
    CSVfile.close()

    # print(Samples[1, :])
    # print(min(Samples[1, :]))
    # print(max(Samples[1, :]))
    # plt.hist(Samples[1, :])
    # plt.show()

# srun --mem=50G --time=00:15:00 python RootToCSV.py
