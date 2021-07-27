import matplotlib.pyplot as plt
import numpy as np
import uproot


def readMultipleRoot(fileName):
    print("Read in:", fileName)
    file = uproot.open(fileName)

    # print(file.classnames())

    tree = file["Detector Data 0"]
    # tree.show()

    Branches = tree.keys()
    NumBranches = len(Branches)

    Temp = tree["Gun_energy_MeV"].array(library="np")
    NumPoints = len(Temp)
    print("Number of Data points", "{:,}".format(NumPoints))

    Data = np.zeros((NumBranches-1, NumPoints))

    Data[0] = Temp
    print("Data[0] = Gun_energy_MeV")

    for i in range(NumBranches-2):
        print("Data[" + str(i+1) + "] = " + Branches[i])
        Data[i+1] = tree[Branches[i]].array(library="np")

    return Data  # 0: Gun_energy; 1: SiVol_0_Edep; 2: SiVol_0_Esec; 3: SiVol_1_Edep; 4: SiVol_1_Esec etc...


if __name__ == "__main__":
    DataTest = readMultipleRoot("/home/anton/Desktop/triton_work/3D/SiChipTestModular/root/1mmal3dmultichip2e9electrons500kev.root")
    print(type(DataTest))
    print(type(DataTest[0][0]))
    plt.hist(DataTest[0], weights=DataTest[1], bins=100)
    plt.yscale("log")
    plt.xlabel("Primary particle energy [MeV]")
    plt.ylabel("Dose per primary energy bin")
    plt.show()
    # plt.savefig("out/MulasProton2e9AlFull-16mm.eps", format='eps')

# srun --mem-per-cpu=50G --time=00:15:00 python ReadG4root.py
