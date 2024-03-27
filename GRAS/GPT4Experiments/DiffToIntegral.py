from GRAS.Read.ReadSpenvis_tri import readSpenvis_tri
import matplotlib.pyplot as plt
import numpy as np

## Read in GTO spectra
GTO_file = "/l/triton_work/Spectra/GTO/spenvis_tri.txt"
Protons, Electrons = readSpenvis_tri(GTO_file)


# Calculate Integral Flux values based on the Differential Flux values
NumberParticlesBetweenTwoEnergyValues = []
for i in range(len(Protons['Energy'])-1):

    # Calulate the number of particles between two energy values
    E1 = Protons['Energy'][i]
    E2 = Protons['Energy'][i+1]

    Diff1 = Protons['Differential'][i]
    Diff2 = Protons['Differential'][i+1]

    # Linear interpolation between the two energy values
    # A = 0.5 * (Diff1 + Diff2) * (E2 - E1)
    # NumberParticlesBetweenTwoEnergyValues.append(A)

    # Alternatively power law interpolation
    # Calculate power-law parameters
    k = np.log(Diff2/Diff1) / np.log(E2/E1)
    a = Diff1 / (E1**k)
    
    # Calculate the integral for this segment
    if k != -1:
        J_segment = (a / (k + 1)) * (E2**(k + 1) - E1**(k + 1))
    else:  # Handle the k = -1 case if necessary
        J_segment = a * np.log(E2 / E1)
    
    NumberParticlesBetweenTwoEnergyValues.append(J_segment)

# Calculate the Integral Flux values
IntegralFlux = []
for i in range(len(NumberParticlesBetweenTwoEnergyValues)):
    # The integral flux is the sum of all particles above the energy value 
    IntegralFlux.append(sum(NumberParticlesBetweenTwoEnergyValues[i:]))

IntegralFlux.append(0)  # Add the last value to the list


# Print the Integral Flux values and relative difference
for i in range(len(Protons['Energy'])):
    print(f"{Protons['Energy'][i]:.4f}", f"{Protons['Integral'][i]:.4e}", f"{IntegralFlux[i]:.4e}", f"{(Protons['Integral'][i] - IntegralFlux[i])/Protons['Integral'][i]:.4f}")

# Plot the Integral Flux values
plt.figure(0)
plt.plot(Protons['Energy'], Protons['Integral'], label='Given Integral Flux')
plt.plot(Protons['Energy'], IntegralFlux, label='Calculated Integral Flux from Differential Flux')
plt.title('Calculated Integral Flux for GTO Protons')
plt.xscale("log")
plt.yscale("log")
plt.xlabel('Energy (MeV)')
plt.ylabel('J(E)')
plt.grid(which="both")
plt.legend()
plt.show()
