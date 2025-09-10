import matplotlib.pyplot as plt

# Electron data
ElectronEnergy = [0.001, 0.01, 0.1, 1, 2.5, 5, 10, 20, 40, 80, 160, 320, 640, 1000, 10000]
ElectronLET = [2000, 401.81, 77.47, 35.79, 37.293, 40.854, 46.85, 57.539, 78.11, 119.26, 202.57, 370.71, 708.86, 1090.7, 10874]


# Proton data
ProtonEnergy = [0.01, 0.1, 1, 2, 5, 10, 25, 50, 100, 150, 200, 500, 1000, 10000, 100000, 1e6, 1e7]
ProtonLET = [7786.6, 11782, 4088, 2602.7, 1359.4, 805.05, 392.61, 229.38, 136, 102.14, 84.543, 52.008, 42.068, 42.592, 51.965, 64.612, 137.44]

plt.figure(figsize=[10,8])

plt.loglog(ElectronEnergy, ElectronLET, label='Electron')
plt.loglog(ProtonEnergy, ProtonLET, label='Proton')

plt.title('LET vs Energy for Electrons and Protons')
plt.xlabel('Energy (MeV)')
plt.ylabel('LET (MeV/cm)')

plt.legend()
plt.grid(True, which="both", ls="--")
#plt.show()
plt.savefig("/l/triton_work/Chess1 GNSS SEE Analysis/Stopping Powers and Ranges/Geant4 Stopping Power and Ranges.pdf", format='pdf', bbox_inches="tight")