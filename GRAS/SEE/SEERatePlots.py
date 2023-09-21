import numpy as np
from IPython.display import display
import matplotlib.pyplot as plt
import pandas as pd

plt.rc('axes', axisbelow=True)

df = pd.read_csv('/l/triton_work/CARRINGTON/SEERates.csv')
FitParams = "Uncorrectable"

DataNames = [#"AP8", "SolarP", "AE8", "CosmicP", "CosmicFe"]
    "Carrington SEP +2 Sigma",
    "Carrington SEP EVT",
    "Carrington SEP -2 Sigma",
    "2003 SPE",
    #"AE9 GTO trapped Electrons",
    "AP9 GTO trapped Protons",
    #"AE9 LEO trapped Electrons",
    "AP9 LEO trapped Protons"]

Colours = ['C1', 'C0', 'C2', 'C8', 'C3', 'C9']

#print(df.to_string())

print(df.keys())

if FitParams == "Hard":
    CrossectionName = "Hard reset"
    C = df[df['Crossection'].str.contains('Hard')]
elif FitParams == "Soft":
    CrossectionName = "Soft reset"
    C = df[df['Crossection'].str.contains('Soft')]
elif FitParams == "Correctable":
    CrossectionName = "LSRAM correctable SBU"
    C = df[~df['Crossection'].str.contains('un')]  # Correctable
elif FitParams == "Uncorrectable":
    CrossectionName = "LSRAM uncorrectable SBU"
    C = df[df['Crossection'].str.contains('un')]

# Initialize dictionaries to hold sums and error sums
#sums = {}
#sum_errors = {}


plt.fill_between(C.Shielding[C.Data == 'Carrington SEP EVT'], C.SEE_Rate[C.Data == 'Carrington SEP EVT'], C.SEE_Rate[C.Data == 'Carrington SEP +2 Sigma'], color='C1', alpha=0.5)
plt.fill_between(C.Shielding[C.Data == 'Carrington SEP EVT'], C.SEE_Rate[C.Data == 'Carrington SEP EVT'], C.SEE_Rate[C.Data == 'Carrington SEP -2 Sigma'], color='C2', alpha=0.5)


i = 0
for DataName in DataNames:
    Cd = C[C.Data == DataName]
    plt.errorbar(Cd.Shielding, Cd.SEE_Rate, yerr=Cd.SEE_Error, capsize=5, elinewidth=1, capthick=2, label=DataName, color=Colours[i])

    i += 1

#plt.plot([1/8e+3]*5, '--', label="1 Bitflip per kB per Second", color='black')
#plt.plot([1/8e+6]*5, ':', label="1 Bitflip per MByte per second", color='black')
plt.plot([1/8e+9]*5, '-.', label="1 Bitflip per GByte per second", color='black')
plt.plot([1/8e+12]*5, '--', label="1 Bitflip per TByte per second", color='black')

#plt.ylim(1e-15, 1e-3)
plt.yscale("log")
plt.grid()
plt.title(CrossectionName)

plt.xlabel("Aluminium Shielding Thickness")
plt.ylabel("Single Bit Upset Rate [s-1 bit-1]")
plt.legend()

plt.savefig("/l/triton_work/CARRINGTON/" + CrossectionName + "Rates.pdf", format='pdf', bbox_inches="tight")