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
    print("Branches:", Branches)
    NumBranches = len(Branches)

    GunEnergy = tree["Gun_energy_MeV"].array(library="np")
    NumPoints = len(GunEnergy)
    print("Number of Data points", "{:,}".format(NumPoints))

    Data = np.zeros((NumBranches, NumPoints))

    for i in range(NumBranches):
        Data[i] = tree[Branches[i]].array(library="np")

    return Data  # ['Sivol_0_Edep_MeV', ... , 'Sivol_0_Esec_MeV', ... , 'Gun_energy_MeV', 'Gun_angle_deg']


if __name__ == "__main__":
    DataTest = readMultipleRoot("/home/anton/Desktop/triton_work/3D/MultiChipTest/root/1-mmal3dmultichip2e9electrons500kev.root")
    print(type(DataTest))
    print(type(DataTest[0]))
    plt.hist(DataTest[-2], weights=DataTest[1], bins=100)
    plt.yscale("log")
    plt.xlabel("Primary particle energy [MeV]")
    plt.ylabel("Dose per primary energy bin")
    plt.show()
    # plt.savefig("out/MulasProton2e9AlFull-16mm.eps", format='eps')

# srun --mem-per-cpu=50G --time=00:15:00 python ReadG4root.py
