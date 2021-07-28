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
    NumVols = int((len(Branches)/2)-1)

    GunEnergy = tree["Gun_energy_MeV"].array(library="np")
    NumPoints = len(GunEnergy)
    print("Number of Data points", "{:,}".format(NumPoints))

    Edep = np.zeros((NumVols, NumPoints))
    Esec = np.zeros((NumVols, NumPoints))

    for i in range(NumVols):
        Edep[i] = tree[Branches[i]].array(library="np")
        Esec[i] = tree[Branches[i+NumVols]].array(library="np")

    return GunEnergy, Edep, Esec  # 0: Gun_energy; 1: SiVol_0_Edep; 2: SiVol_0_Esec; 3: SiVol_1_Edep; 4: SiVol_1_Esec etc...


if __name__ == "__main__":
    GunEnergyTest, EdepTest, EsecTest = readMultipleRoot("/home/anton/Desktop/triton_work/3D/SiChipTestModular/root/1mmal3dmultichip2e9electrons500kev.root")
    print(type(GunEnergyTest))
    print(type(GunEnergyTest[0]))
    plt.hist(GunEnergyTest, weights=EdepTest[0], bins=100)
    plt.yscale("log")
    plt.xlabel("Primary particle energy [MeV]")
    plt.ylabel("Dose per primary energy bin")
    plt.show()
    # plt.savefig("out/MulasProton2e9AlFull-16mm.eps", format='eps')

# srun --mem-per-cpu=50G --time=00:15:00 python ReadG4root.py
