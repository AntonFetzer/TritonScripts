import numpy as np
import matplotlib.pyplot as plt
from GRAS.Read.ReadSD2Q import readSDQ2
from GRAS.Dependencies.TotalDose import totalDose
from GRAS.Dependencies.MergeTotalDose import mergeTotalDose


##################################### Reading in GRAS data ########################################

Path = "/l/triton_work/Shielding_Curves/FS1-SunSensors/"

#### FS1 ######
Electrons = totalDose(Path + "AE9-FS1-mission/Res")
Protons = totalDose(Path + "AP9-FS1-mission/Res")

SolarProtons = totalDose(Path + "SAPPHIRE-FS1/Res")

Total = mergeTotalDose([Electrons, Protons, SolarProtons])


############################# Plotting GRAS data ########################################
plt.figure(0, [5, 5])
# Aluminium Shielding Thickness [mm]
# 101 Tiles from 0 to 10 mm
x = np.linspace(0, 10, num=101, dtype=float, endpoint=True)
for i in x:
    print(i)

# Plot FS1 total
plt.errorbar(x, Total['dose'], yerr=Total['error'], capsize=5, linestyle='', color='C8', label="GRAS Total Dose")
# Plot FS1 Electrons
plt.errorbar(x, Electrons['dose'], yerr=Electrons['error'], capsize=5, linestyle='', color='C0', label="GRAS Electrons")
# Plot FS1 Protons
plt.errorbar(x, Protons['dose'], yerr=Protons['error'], capsize=5, linestyle='', color='C1', label="GRAS Trapped Protons")
# Plot FS1 Solar Protons
plt.errorbar(x, SolarProtons['dose'], yerr=SolarProtons['error'], capsize=5, linestyle='', color='C2', label="GRAS Solar Protons")



####### Plot 1kRad line #########
CriticalDose = [1 for i in x]
plt.plot(x, CriticalDose, color='k', linewidth=2, label='1 krad')
CriticalDose = [10 for i in x]
plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='10 krad')
# CriticalDose = [100 for i in x]
# plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='100 krad')


##### SHIELDOSE #####
# ------------------------------- Import and Plot SHIELDOSE Data -------------------------------------------------------
SDData = readSDQ2("/l/triton_work/Spectra/FS1/SHIELDOSE-2Q/spenvis_sqoCustomThicknesses.txt")
# SDDataCollumns = ['Thickness', 'Total','Electrons','Bremsstrahlung','Trapped Protons','Solar Protons']
plt.plot(SDData['Thickness'], SDData['Total'], 'C8.', label="SHIELDOSE-2Q Total Dose")
plt.plot(SDData['Thickness'], SDData['Electrons'] + SDData['Bremsstrahlung'], 'C0.', label="SHIELDOSE-2Q Electrons")
plt.plot(SDData['Thickness'], SDData['Trapped Protons'], 'C1.', label="SHIELDOSE-2Q Trapped Protons")
plt.plot(SDData['Thickness'], SDData['Solar Protons'], 'C2.', label="SHIELDOSE-2Q Solar Protons")



plt.title("Ionising Dose Behind Shielding")
plt.xlabel("Aluminium Shielding Thickness [mm]")
plt.ylabel("Total Ionising Dose [krad]")
plt.yscale("log")
# plt.xscale("log")


plt.xlim(-0.5, 5)
# plt.ylim(2e-1, 2e+2)
plt.grid(which='both')


################### Re-ordering the legend labels ############################
# Store the handles and labels when plotting
# handles, labels = plt.gca().get_legend_handles_labels()

# # Reorder the handles and labels to your desired order
# # Example: Suppose you want to reorder them based on some criteria or the original plotting order
# order = [0, 1, 4, 5, 2, 6, 3]  # specify the desired order of indices

# # Apply the reordering
# handles = [handles[i] for i in order]
# labels = [labels[i] for i in order]

# # Then pass them to plt.legend()
# plt.legend(handles, labels)
plt.legend()


plt.savefig(Path + "CompareGRASShieldose.pdf", format='pdf', bbox_inches="tight")
# plt.show()

# Print table comparing SHIELDOSE and GRAS doses to CSV file
CSVFile = open("/l/triton_work/Shielding_Curves/FS1-SunSensors/CompareGRASShieldose.csv", 'w')

CSVFile.write("x, GRAS Total, GRAS Electrons, GRAS Trapped Protons, GRAS Solar Protons, SHIELDOSE Thickness, SHIELDOSE Total, SHIELDOSE Electrons, SHIELDOSE Trapped Protons, SHIELDOSE Solar Protons\n")

for i in range(len(x)):
    CSVFile.write(f"{x[i+1]}, {Total['dose'][i+1]}, {Electrons['dose'][i+1]}, {Protons['dose'][i+1]}, {SolarProtons['dose'][i+1]}, {SDData['Thickness'][i]}, {SDData['Total'][i]}, {SDData['Electrons'][i] + SDData['Bremsstrahlung'][i]}, {SDData['Trapped Protons'][i]}, {SDData['Solar Protons'][i]}\n")