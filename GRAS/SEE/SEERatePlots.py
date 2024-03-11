import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rc('axes', axisbelow=True)

df = pd.read_csv('/l/triton_work/LET/Foresail1-Hercules/SEERates.csv')
FitParams = "SEFI"

DataNames = [
    "FS1-Solar Protons",
    "FS1-Cosmic Protons",
    "FS1-Cosmic Helium",
    "FS1-Cosmic Iron",
    "FS1-Cosmic Oxygen",
    "FS1-Trapped Protons",
    "FS1-Trapped Electrons"
    ]

Colours = ['C1', 'C0', 'C2', 'C8', 'C3', 'C9', 'C4']

#print(df.to_string())

print(df.keys())

if FitParams == "SEU":
    CrossectionName = "Hercules SEU"
    C = df[df['Crossection'].str.contains('SEU')]
elif FitParams ==  "SET":
    CrossectionName = "Hercules SET"
    C = df[df['Crossection'].str.contains('SET')]
elif FitParams == "SEFI":
    CrossectionName = "Hercules SEFI"
    C = df[df['Crossection'].str.contains('SEFI')]


for i, DataName in enumerate(DataNames):
    Cd = C[C.Data == DataName]
    plt.errorbar(Cd.Shielding, Cd.SEE_Rate, yerr=Cd.SEE_Error, capsize=5, elinewidth=1, capthick=2, label=DataName, color=Colours[i])

#plt.plot([1/8e+3]*5, '--', label="1 Bitflip per kB per Second", color='black')
#plt.plot([1/8e+6]*5, ':', label="1 Bitflip per MByte per second", color='black')
#plt.plot([1/8e+9]*5, '-.', label="1 Bitflip per GByte per second", color='black')
#plt.plot([1/8e+12]*5, '--', label="1 Bitflip per TByte per second", color='black')

#plt.ylim(1e-15, 1e-3)
plt.yscale("log")
plt.grid()
plt.title(CrossectionName)

plt.xlabel("Aluminium Shielding Thickness")
plt.ylabel("Single Bit Upset Rate [s-1 bit-1]")
plt.legend(loc='lower left')

#plt.savefig("/l/triton_work/CARRINGTON/" + CrossectionName + "Rates.pdf", format='pdf', bbox_inches="tight")
plt.show()