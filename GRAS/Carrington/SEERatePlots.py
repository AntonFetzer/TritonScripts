import numpy as np
from IPython.display import display
import matplotlib.pyplot as plt
import pandas as pd

plt.rc('axes', axisbelow=True)

df = pd.read_csv('/l/TritonPlots/Luna/SEE/SEERates.csv')
CorrectableOrNot = 0  # 0 for correctable 1 for uncorrectable

DataNames = [  # "Carrington SEP +2 Sigma",
    # "Carrington SEP EVT",
    # "Carrington SEP -2 Sigma",
    # "2003 SPE",
    "AP9 GTO trapped protons",
    "AP9 LEO trapped protons",
    "Cosmic Protons",
    "Solar Protons"
]

Colours = ['C1', 'C0', 'C2', 'C8', 'C3', 'C9']

#print(df.to_string())

print(df.keys())

if CorrectableOrNot == 0:
    C = df[~df['Crossection'].str.contains('un')]  # Correctable
elif CorrectableOrNot == 1:
    C = df[df['Crossection'].str.contains('un')]

#plt.fill_between(C.Shielding[C.Data == 'Carrington SEP EVT'], C.SEE_Rate[C.Data == 'Carrington SEP EVT'], C.SEE_Rate[C.Data == 'Carrington SEP +2 Sigma'], color='C1', alpha=0.5)
#plt.fill_between(C.Shielding[C.Data == 'Carrington SEP EVT'], C.SEE_Rate[C.Data == 'Carrington SEP EVT'], C.SEE_Rate[C.Data == 'Carrington SEP -2 Sigma'], color='C2', alpha=0.5)

i=0
for DataName in DataNames:
    Cd = C[C.Data == DataName]
    plt.errorbar(Cd.Shielding, Cd.SEE_Rate, yerr=Cd.SEE_Error, capsize=5, elinewidth=1, capthick=2, label=DataName, color=Colours[i])
    i += 1

plt.fill_between(Cd.Shielding[C.Data == 'Carrington SEP EVT'], Cd.SEE_Rate[C.Data == 'Carrington SEP EVT'], Cd.SEE_Rate[C.Data == 'Carrington SEP +2 Sigma'], color='C1', alpha=0.5)
#plt.fill_between(Cd.Shielding, Cd.SEE_Rate,, color='C2', alpha=0.5)

#plt.plot([1/8e+3]*6, '--', label="1 Bitflip per kB per Second", color='black')
plt.plot([1/8e+6]*6, ':', label="1 Bitflip per MByte per second", color='black')
plt.plot([1/8e+9]*6, '-.', label="1 Bitflip per GByte per second", color='black')
#plt.plot([1/8e+12]*6, '--', label="1 Bitflip per TByte per second", color='black')

#plt.ylim(1e-15, 1e-3)
plt.yscale("log")
plt.grid()
if CorrectableOrNot == 0:
    plt.title("Correctable Single Bit Upset Rates")  # Correctable
elif CorrectableOrNot == 1:
    plt.title("Uncorrectable Single Bit Upset Rates")

plt.xlabel("Aluminium Shielding Thickness")
plt.ylabel("Single Bit Upset Rate [s-1 bit-1]")
plt.legend()

if CorrectableOrNot == 0:
    plt.savefig("/l/TritonPlots/Luna/CorrectableSEERate.svg", format='svg', bbox_inches="tight")  # Correctable
elif CorrectableOrNot == 1:
    plt.savefig("/l/TritonPlots/Luna/UncorrectableSEERate.svg", format='svg', bbox_inches="tight")

