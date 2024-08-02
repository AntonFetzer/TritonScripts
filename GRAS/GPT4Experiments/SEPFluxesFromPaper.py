import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammaincc


# Parameters for Hydrogen on 10/28/03 for the Ellison-Ramaty formula
K = 1.35e9  # Normalization
gamma = -1.05  # Gamma
E0 = 28.2  # E0 in MeV/nuc

# Implementing the Ellison-Ramaty spectral shape function for differential flux
def ERdiff(E, K, gamma, E0):
    return K * E**gamma * np.exp(-E/E0)

# Creating an array of energies for plotting
E_values = np.geomspace(1e-1, 600, 100)  # from 0.1 to 600 MeV

# Calculating the differential flux values for the given energies
ERdiffValues = np.array([ERdiff(E, K, gamma, E0) for E in E_values])

#Print Energies and Differential Flux values
for i in range(len(E_values)):
    print(f"{E_values[i]:.4f}", f"{ERdiffValues[i]:.4f}")

# Plotting the differential hydrogen flux data for 10/28/03
plt.figure(0)
plt.loglog(E_values, ERdiffValues, label='Ellison-Ramaty 10/28/03')
plt.title('Ellison-Ramaty Spectra for Hydrogen on Multiple Events')
plt.xlabel('Energy (MeV/nuc)')
plt.ylabel('dJ/dE')
plt.grid(which="both")
plt.ylim(1e0, 1e10)
plt.legend()


def ERIntegral(E, K, gamma, E0):
    s = 1 - gamma  # s parameter for the incomplete gamma function
    x = E / E0     # x parameter for the incomplete gamma function
    # The upper incomplete gamma function gammaincc takes arguments (s, x)
    # and returns the normalized upper incomplete gamma function, which we need to multiply by gamma(s)
    return K * E0**(1 - gamma) * gammaincc(s, x) * np.math.gamma(s)

# Calculating the integral flux values for the given energies
integral_flux_values_explicit = [ERIntegral(E, K, gamma, E0) for E in E_values]

# Plot the integral flux computed using the explicit formula

plt.figure(1)
plt.loglog(E_values, integral_flux_values_explicit, label='Integral Flux (Explicit) 10/28/03', linestyle='--')
plt.title('Integral Flux for Hydrogen on 10/28/03 (Explicit Formula)')
plt.xlabel('Energy (MeV/nuc)')
plt.ylabel('Integral Flux')
plt.grid(which="both")
plt.ylim(1e0, 1e10)
plt.legend()

plt.show()
