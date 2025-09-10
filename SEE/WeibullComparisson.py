import numpy as np
import matplotlib.pyplot as plt

WD = {}

WD["PolarFireCorrectable"] = {
    "L0": 0.4, 
    "W": 18,
    "S": 0.98,
    "A0": 4.01e-9,
    "color": 'C0'
}
WD["DUT6"] = {
    "L0": 2.801, 
    "W": 37.833,
    "S": 2.972,
    "A0": 5.185E-09,
    "color": 'C1'
}
WD["DUT7"] = {
    "L0": 0.851, 
    "W": 35.460,
    "S": 5.575,
    "A0": 5.205E-09,
    "color": 'C2'
}
# W["DUT6+7"] = {
#     "L0": 0.11214, 
#     "W": 36.4286,
#     "S": 4.44737,
#     "A0": 5.1852E-09
# }
# W["nanoXplore "] = {
#     "L0": 0.11, 
#     "W": 36,
#     "S": 4.4,
#     "A0": 5.2E-09
# }
WD["PolarFireUncorrectable"] = {
    "L0": 0.4,
    "W": 1,
    "S": 0.4,
    "A0": 5.5e-14,
    "color": 'C3'
}


### Data points
WD["PolarFireCorrectable"]["LET"] = [1.285140562248996, 2.831325301204819, 5.291164658634538, 9.29718875502008, 23.84538152610442]
WD["PolarFireCorrectable"]["Rate"] = [1.5262399906547595e-10, 7.644950712394795e-10, 7.911576418029598e-10, 7.777121046581536e-10, 3.012369131206531e-9]

WD["PolarFireUncorrectable"]["LET"] = [1.4257028112449799, 2.831325301204819, 5.431726907630522, 9.29718875502008, 23.84538152610442]
WD["PolarFireUncorrectable"]["Rate"] = [4.6150786268561354e-14, 4.6150786268561354e-14, 4.536453022869048e-14, 4.695066966348511e-14, 4.536453022869048e-14]

WD["DUT6"]["LET"] = [3.3, 10, 20.4, 32.4, 62.5, 3.3, 20.4, 32.4, 62.5, 3.3, 10, 20.4, 32.4]
WD["DUT6"]["Rate"] = [3E-14, 8E-12, 6E-10, 2E-09, 5E-09, 2E-13, 2E-09, 4E-09, 1E-08, 1E-11, 2E-09, 3E-09, 5E-09]

WD["DUT7"]["LET"] = [5.7, 16, 32.4, 62.5, 5.7, 32.4, 62.5]
WD["DUT7"]["Rate"] = [4E-14, 3E-10, 2E-09, 5E-09, 3E-10, 5E-09, 2E-08]

# def weibull_cdf(x, L0, W, S, A0):
    # return 1 - np.exp(-((x - L0) / W) ** S * A0)

    
def f(LET, L0, W, S, A0):
    """
    Parameters:
    LET (numpy.ndarray): Array of LET values.

    Returns:
    numpy.ndarray: Array of rate estimates.
    """
    
    result = np.zeros_like(LET)
    # Remove all values of LET that are less than x0
    mask = LET > L0
    # A * (1 - np.exp(-((x - x0) / W) ** s))
    result[mask] = A0 * (1 - np.exp(-((LET[mask] - L0) / W) ** S))

    return result

x = np.geomspace(1e-1, 100, 1000)
#print(x)

plt.figure()

# WD.pop("PolarFireCorrectable")
# WD.pop("PolarFireUncorrectable")
# WD.pop("DUT6")
# WD.pop("DUT7")

for key in WD:
    print(key)
    L0 = WD[key]["L0"]
    W = WD[key]["W"]
    S = WD[key]["S"]
    A0 = WD[key]["A0"]
    # y = weibull_cdf(x, L0, W, S, A0)
    y = f(x, L0, W, S, A0)
    plt.plot(x, y, label=key, color=WD[key]["color"])

# Plot data points
for key in WD:
    plt.scatter(WD[key]["LET"], WD[key]["Rate"], color=WD[key]["color"])

plt.xlabel('LET [MeV cm$^2$/mg]')
plt.ylabel('Cross section cm$^{-2}$')
plt.title('Weibull CDF Comparison')
plt.legend(loc='upper left')
plt.grid(True)
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-1, 80)
plt.ylim(1e-14, 1e-8)
#plt.show()
plt.savefig("/l/triton_work/Python/GRAS/SEE/WeibullCDFComparison.pdf", format='pdf', bbox_inches='tight')