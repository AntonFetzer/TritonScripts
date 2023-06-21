import numpy as np
from IPython.display import display
import matplotlib.pyplot as plt
import pandas as pd

plt.rc('axes', axisbelow=True)

df = pd.read_csv('/home/anton/Desktop/triton_work/Chess1 GNSS SEE Analysis/SEERates.csv')
HardSoft = 1  # 0 for Soft 1 for Hard

DataNames = ["AP8", "SolarP", "AE8", "CosmicP", "CosmicFe"]

Colours = ['C1', 'C0', 'C2', 'C8', 'C3', 'C9']

#print(df.to_string())

print(df.keys())

if HardSoft == 0:
    C = df[df['Crossection'].str.contains('Soft')]  # Correctable
elif HardSoft == 1:
    C = df[df['Crossection'].str.contains('Hard')]


# Initialize dictionaries to hold sums and error sums
sums = {}
sum_errors = {}

i = 0
for DataName in DataNames:
    Cd = C[C.Data == DataName]
    plt.errorbar(Cd.Shielding, Cd.SEE_Rate, yerr=Cd.SEE_Error, capsize=5, elinewidth=1, capthick=2, label=DataName, color=Colours[i])

    # Add up SEE_Rate values and squared SEE_Error values for each shielding thickness
    for index, row in Cd.iterrows():
        if row['Shielding'] in sums:
            sums[row['Shielding']] += row['SEE_Rate']
            sum_errors[row['Shielding']] += row['SEE_Error']**2
        else:
            sums[row['Shielding']] = row['SEE_Rate']
            sum_errors[row['Shielding']] = row['SEE_Error']**2
    i += 1

# Convert sums dictionary to a DataFrame and plot it
sums_df = pd.DataFrame(list(sums.items()), columns=['Shielding', 'SEE_Rate'])
sum_errors_df = pd.DataFrame(list(sum_errors.items()), columns=['Shielding', 'SEE_Error'])

# Take square root of sum_errors to get total error
sum_errors_df['SEE_Error'] = np.sqrt(sum_errors_df['SEE_Error'])

# Plot total with error bars
plt.errorbar(sums_df.Shielding, sums_df.SEE_Rate, yerr=sum_errors_df['SEE_Error'], capsize=5, elinewidth=1, capthick=2, label='Total', color='C4')

#plt.plot([1/8e+3]*6, '--', label="1 Bitflip per kB per Second", color='black')
#plt.plot([1/8e+6]*6, ':', label="1 Bitflip per MByte per second", color='black')
#plt.plot([1/8e+9]*6, '-.', label="1 Bitflip per GByte per second", color='black')
#plt.plot([1/8e+12]*6, '--', label="1 Bitflip per TByte per second", color='black')

#plt.ylim(1e-15, 1e-3)
plt.yscale("log")
plt.grid()
if HardSoft == 0:
    plt.title("Soft Errors")  # Correctable
elif HardSoft == 1:
    plt.title("Hard Errors")

plt.xlabel("Aluminium Shielding Thickness")
plt.ylabel("Total Error Number")
plt.legend()

if HardSoft == 0:
    plt.savefig("/home/anton/Desktop/triton_work/Chess1 GNSS SEE Analysis/Soft Error Rates.pdf", format='pdf', bbox_inches="tight")  # Correctable
elif HardSoft == 1:
    plt.savefig("/home/anton/Desktop/triton_work/Chess1 GNSS SEE Analysis/Hard Error Rates.pdf", format='pdf', bbox_inches="tight")