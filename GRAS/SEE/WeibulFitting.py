import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# The original Weibull function given
def weibull_func(x, A, x0, W):
    return A * (1 - np.exp(-((x - x0) / W)))

HardSoft = 1  # 0 for Soft 1 for Hard
# Your data points
x_data = np.array([0.00984463519313305, 0.00583690987124464, 0.00438369098712446, 0.00362845493562232])   # LET [MeV cm2 mg-1]
if HardSoft == 0:
    y_data = np.array([2.37388724035608E-09, 1.84119677790564E-09, 2.14564369310793E-09, 2.68987341772152E-09])   # Soft reset cross section [cm-2]
elif HardSoft == 1:
    y_data = np.array([1.00148367952522E-09, 3.73993095512083E-10, 7.47724317295189E-10, 7.51582278481013E-10])   # Hard reset cross section [cm-2]

p0 = [np.max(y_data), 0, np.std(x_data)]
bounds = ([0, 0, 0], [np.inf, np.inf, np.inf])  # set different bounds for each parameter
# Fit the Weibull function to the data
popt, pcov = curve_fit(weibull_func, x_data, y_data, p0=p0, bounds=bounds,  maxfev=5000)

# Print the optimal parameters
A, x0, W = popt
print(f"A = {A}, x0 = {x0}, W = {W}")

# Plot the data and the fitted function
plt.scatter(x_data, y_data, label='Data')
plt.plot(np.linspace(0, max(x_data)*1.5, 500), weibull_func(np.linspace(0, max(x_data)*1.5, 500), *popt),
         'r-', label=f'Fit: A={A:.3e}, x0={x0:.3e}, W={W:.3e}')

plt.xlabel('LET [MeV cm2 mg-1]')
if HardSoft == 0:
    plt.ylabel('# Total Soft Reset cross section [cm-2] ')
    plt.title("Soft Errors cross section U-Blox ZED-F9P GNSS Receiver at PSI")
elif HardSoft == 1:
    plt.ylabel('# Total Hard Reset cross section [cm-2] ')
    plt.title("Hard Errors cross section U-Blox ZED-F9P GNSS Receiver at PSI")
plt.xlim(0, max(x_data)*1.5)
plt.ylim(0, max(y_data)*1.5)
plt.grid('both')
plt.legend()

if HardSoft == 0:
    plt.savefig("/l/triton_work/Chess1 GNSS SEE Analysis/Fits/Soft Resets.pdf", format='pdf', bbox_inches="tight")
elif HardSoft == 1:
    plt.savefig("/l/triton_work/Chess1 GNSS SEE Analysis/Fits/Hard Resets.pdf", format='pdf', bbox_inches="tight")


#plt.show()
