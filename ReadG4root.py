import matplotlib.pyplot as plt
import numpy as np
import uproot


# ---------------------- Dont Modify -------- Dont break compatability -----------------------

# ---------------------- Dont Modify -------- Dont break compatability -----------------------

# ---------------------- Dont Modify -------- Dont break compatability -----------------------

def readG4root(fileName):
    print("Read in:", fileName)
    file = uproot.open(fileName)

    # print(file.classnames())

    tree = file["Detector Data 0"]
    # tree.show()

    Temp = tree["Gun_energy_MeV"].array(library="np")

    Data = np.zeros((3, len(Temp)))

    Data[0] = Temp
    Data[1] = tree["Sivol_Edep_MeV"].array(library="np")
    Data[2] = tree["Sivol_Esec_MeV"].array(library="np")

    return Data  # 0: Gun_energy; 1: SiVol_Edep; 2: SiVol_Esec


if __name__ == "__main__":
    DataTest = readG4root(
        "/home/anton/Desktop/triton_work/TEST/SiChip in 1cm cube for Spenvis Test/root/electrons500kev1mm2SiChip.root")
    print(type(DataTest))
    print(type(DataTest[0][0]))
    plt.hist(DataTest[0], bins=1000)
    plt.yscale("log")
    plt.xlabel("Primary particle energy [MeV]")
    plt.ylabel("Number of counts per bin")
    plt.show()
    # plt.savefig("out/MulasProton2e9AlFull-16mm.eps", format='eps')

# srun --mem-per-cpu=50G --time=00:15:00 python ReadG4root.py
