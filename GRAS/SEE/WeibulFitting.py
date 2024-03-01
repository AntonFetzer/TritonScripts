import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Weibull function
# A: amplitude
# x0: threshold
# W: width
# s: shape
def weibull_func(x, A, x0, W, s):
    return A * (1 - np.exp(-((x - x0) / W) ** s))

# Your data points
# x_data = np.array([30.260, 3.593, 0.965])   # LET [MeV cm2 mg-1]
# y_data = np.array([7.505e-12, 1.138e-11, 1.418e-12])   # cross section [cm-2]

# x_data = np.array([3.576, 30.243, 48.165])   # LET [MeV cm2 mg-1]
# y_data = np.array([3.656E-07, 2.548E-07, 2.269E-07])   # cross section [cm-2]

x_data = np.array([0.963, 2.246, 3.479, 6.123])   # LET [MeV cm2 mg-1]
y_data = np.array([9.031E-08, 2.632E-07, 1.447E-06, 2.254E-06])   # cross section [cm-2]

# Initial guess for the parameters
# p0 = [A, x0, W, s]
# A = max(y_data) because the maximum value of the Weibull function is close to A
# x0 = min(x_data)*99 because the threshold is close to but below the minimum value of the x_data
# W = std(x_data) because the width is on the order of the standard deviation of the x_data ?
# s = 1 because 1 is the most stable exponent.
p0 = [np.max(y_data)*1.1, np.min(x_data)*0.99, np.std(x_data), 1]

# Bounds for the parameters
# All parameters are positive
bounds = ([0, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf])  # set different bounds for each parameter

# Fit the Weibull function to the data
# popt: optimal parameters
# pcov: covariance matrix
# maxfev: maximum number of function evaluations
popt, pcov = curve_fit(weibull_func, x_data, y_data, p0=p0, bounds=bounds,  maxfev=5000)

# Print the optimal parameters
A, x0, W, s = popt
print(f"A = {A}, x0 = {x0}, W = {W}, s = {s}")

# Plot the data points
plt.plot(x_data, y_data, '.', label='Data')

# Generate x values for the Weibull fit
# 500 x values between 50% of the minimum x value and 150% of the maximum x value
x_Fit = np.geomspace(min(x_data)*0.5, max(x_data)*1.5, 500)              
y_Fit = weibull_func(x_Fit, *popt)
Fit_label = f'Fit: A={A:.3e}, x0={x0:.3e}, W={W:.3e}, s={s:.3f}'
# Plot the Weibull fit
plt.plot(x_Fit, y_Fit, 'r-', label=Fit_label)

plt.xlabel('LET [MeV cm2 mg-1]')
plt.ylabel('Cross Section [cm-2] ')
plt.title("SEE cross section vs LET")
plt.xlim(min(x_data)*0.5, max(x_data)*1.5)
plt.ylim(min(y_data)*0.2, max(y_data)*1.5)
plt.xscale('log')
plt.yscale('log')
plt.grid('both')
plt.legend()

#plt.show()

# Save the plot
plt.savefig('/l/triton_work/LET/Foresail1-Hercules/Fig13_SEFI_CrossSection_vs_LET_WeibullFitting.png', dpi=300, bbox_inches='tight')
