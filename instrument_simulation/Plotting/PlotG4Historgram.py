import matplotlib.pyplot as plt
from instrument_simulation.Dependencies.ReadMultipleRoot import readMultipleRoot

File = "/l/triton_work/3D/MultiChipTest/TestWithoutShield/root/block-2e9protons10mev.root"

Data = readMultipleRoot(File)
# Data[8] = Gun_energy_MeV
# Data[9] = Gun_angle_deg


#for i in range(10):
#    plt.hist(Data[i], bins=100, label=str(i))

plt.hist(Data[8], weights=Data[8], bins=100, label="Gun_energy_MeV")
plt.hist(Data[8], weights=Data[0], bins=100, label="Sivol_0_Edep_MeV")
plt.hist(Data[8], weights=Data[5], bins=100, label="Sivol_0_Esec_MeV")

plt.yscale("log")
plt.grid(which='major')
plt.title("Histogram of primary energies")
plt.xlabel("Primary particle energy [MeV]")
plt.ylabel("Number of counts per primary energy bin")
plt.legend()
plt.show()
