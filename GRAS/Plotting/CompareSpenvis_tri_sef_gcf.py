from GRAS.Read.ReadSpenvis_tri import readSpenvis_tri
from GRAS.Read.ReadSpenvis_sef import readSpenvis_sef
from GRAS.Read.ReadSpenvis_gcf import readSpenvis_gcf
from GRAS.Read.ReadGPSMacro import readGPSMacro
import numpy as np
import matplotlib.pyplot as plt

TRIProt, TRIElec = readSpenvis_tri("/l/triton_work/Spectra/A9-GTO/spenvis_tri.txt")
SEF = readSpenvis_sef("/l/triton_work/Spectra/SAPPHIRE-GTO/spenvis_sef.txt")
GCF = readSpenvis_gcf("/l/triton_work/Spectra/ISO-GTO/spenvis_gcf.txt")
#GPS = readGPSMacro("/l/triton_work/Spectra/SuperGTO/SuperGTOElectrons.mac")

Species = ['H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl',
           'Ar', 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se',
           'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb',
           'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
           'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At',
           'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U ']

fig1 = plt.figure(1)

####################################################### Trapped Particles
for l in range(np.shape(TRIProt)[1]):
    if TRIProt[1, l] == 0 or TRIProt[2, l] == 0:
        TRIProt = np.delete(TRIProt, l, 1)

for l in range(np.shape(TRIElec)[1]):
    if TRIElec[1, l] == 0 or TRIElec[2, l] == 0:
        TRIElec = np.delete(TRIElec, l, 1)

plt.plot(TRIElec[0], TRIElec[2], label="GT0 Electrons")
plt.plot(TRIProt[0], TRIProt[2], label="GT0 Protons")


####################################################### Solar Particles
#for i in range(np.shape(SEF)[0]):
#    FluxMax = SEF[i, 0, 1]
#    if FluxMax > 50:
#        print("Solar", Species[i], FluxMax)
#        plt.plot(SEF[i, :, 0], SEF[i, :, 1], label="Solar " + Species[i] + " Ions")


plt.plot(SEF[0, :, 0], SEF[0, :, 2], label="Solar Protons")
plt.plot(SEF[1, :, 0], SEF[1, :, 2], label="Solar Helium Ions")
plt.plot(SEF[25, :, 0], SEF[25, :, 2], label="Solar Iron Ions")


####################################################### Cosmic Particles ISO
#for i in range(np.shape(GCF)[0]):
#    FluxMax = GCF[i, 0, 1]
#    if FluxMax > 0.1:
#        print("Cosmic", Species[i], FluxMax)
#        plt.plot(GCF[i, :, 0], GCF[i, :, 1], label="Cosmic " + Species[i] + " Ions")

plt.plot(GCF[0, :, 0], GCF[0, :, 2], label="Cosmic Protons", color='C8')
plt.plot(GCF[1, :, 0], GCF[1, :, 2], label="Cosmic Helium Ions", color='C9')
plt.plot(GCF[25, :, 0], GCF[25, :, 2], label="Cosmic Iron Ions", color='C7')


####################################################### GPS FIles
#plt.plot(GPS[0], GPS[1], 'x', label="GPS Differential Electron Flux")


#plt.xlim(left=0.5)
plt.yscale("log")
plt.xscale("log")
plt.grid(which="both")
plt.title("Differential particle flux")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Differential Flux [cm-2 s-1 MeV-1]")
plt.legend()
#plt.show()

plt.savefig("/l/TritonPlots/Paper/FluxComparison.pdf", format='pdf', bbox_inches="tight")
