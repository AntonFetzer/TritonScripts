import numpy as np
import matplotlib.pyplot as plt
from GRAS.Read.ReadSD2Q import readSDQ2
from GRAS.Dependencies.TotalDose import totalDose
from GRAS.Dependencies.MergeTotalDose import mergeTotalDose

##### SHIELDOSE #####

plt.figure(0, [5, 8])
# ------------------------------- Import and Plot SHIELDOSE Data -------------------------------------------------------
ConvertKrad = 1 / 1000  #  1/1000 to convert to krad
# mmtogcm2 = 0.27  # 0.27 to convert mm to g/cm2 of aliminium

SDData = readSDQ2("/l/triton_work/Spectra/Van-Allen-Belt-Probes/Shieldose/spenvis_sqo.txt")
# SDDataCollumns = ['Aluminium Thickness', 'Total Dose', 'Electrons', 'Bremsstrahlung', 'Protons']
plt.plot(SDData[:, 0], ( SDData[:, 2]+SDData[:, 3] ) * ConvertKrad, 'C0--', label="SHIELDOSE-2Q trapped Electrons")
#plt.plot(SDData[:, 0] * mmtogcm2, SDData[:, 3] * ConvertKrad, 'C1-..', label="SHIELDOSE-2Q Bremsstrahlung")
plt.plot(SDData[:, 0], SDData[:, 4] * ConvertKrad, 'C1-', label="SHIELDOSE-2Q trapped Protons")
plt.plot(SDData[:, 0], SDData[:, 1] * ConvertKrad, 'C2:', label="SHIELDOSE-2Q Total trapped particles")

####### GRAS #######

Path = "/l/triton_work/Shielding_Curves/Carrington/"
res_suffix = "/Res/"

Names = [
    # 'Carrington-SEP-Expected-Diff',
    # 'Carrington-SEP-Expected-Int',

    'Carrington-SEP-Expected-Int-With0',
    'Carrington-SEP-Minus2Sigma-Int-With0',
    'Carrington-SEP-Plus2Sigma-Int-With0',

    # 'CarringtonElectronDiffPow',
    # 'CarringtonElectronDiffPowTabelated',
    'CarringtonElectronINTEGRALPowTabelated',

    # 'GEO-AE9-mission',# Similar but slightly less than VAB-AE9-mission, therefore not interesting for shielding
    # 'GEO-AP9-mission', # The trapped proton spectrum on GEO ends at 8MeV, therefore not intetrsting for shielding

    # 'GEO-SolarProton-mission',
    # 'GEO-SolarProton-5minPeakFlux', # Similar Flux as the carrington protons, but onlf for 5 minutes, this makes the TID negligible

    # 'GEO-CosmicProton-mission',  # Cosmics have super low fluxes --> negligible TID
    # 'GEO-CosmicIron-mission',    # Cosmics have super low fluxes --> negligible TID

    'ISS-AE9-mission',
    'ISS-AP9-mission',

    'VAB-AE9-mission',
    'VAB-AP9-mission',
]

""" 
for name in Names:
   path = Path + name + res_suffix
   ShieldingCurves[name] = totalDose(path)
   ShieldingCurves[name]['color'] = 'C' + str(Names.index(name))
   ShieldingCurves[name]['label'] = name


# The carrington fluxes are per second, but we are interested in the event fluence which corresponds to 1.23 days
for name in Names:
    if "Carrington" in name:
        ShieldingCurves[name]['dose'] *= 1.23 * 24 * 60 * 60
        ShieldingCurves[name]['error'] *= 1.23 * 24 * 60 * 60
    # # The peak fluxes are per second, but we are interested in the event fluence which corresponds to 5 minutes
    # if "PeakFlux" in name:
    #     ShieldingCurves[name]['dose'] *= 5 * 60
    #     ShieldingCurves[name]['error'] *= 5 * 60
    # The ISS fluxes are per month, but they are so low that we better show the 1 year dose
    if "ISS" in name: 
        ShieldingCurves[name]['dose'] *= 12
        ShieldingCurves[name]['error'] *= 12
    # # The GEO solar protons are per month, but TID would be low, so we show the 1 year dose
    # if "GEO-SolarProton" in name:
    #     ShieldingCurves[name]['dose'] *= 120
    #     ShieldingCurves[name]['error'] *= 120

ShieldingCurves['Carrington-SEP-Expected-Int-With0']['label'] = 'Carrington Proton Event Fluence'
ShieldingCurves['Carrington-SEP-Minus2Sigma-Int-With0']['label'] = 'Carrington Proton − 2\u03C3'
ShieldingCurves['Carrington-SEP-Plus2Sigma-Int-With0']['label'] = 'Carrington Proton + 2\u03C3'
ShieldingCurves['Carrington-SEP-Expected-Int-With0']['color'] = 'red'

ShieldingCurves['CarringtonElectronINTEGRALPowTabelated']['label'] = 'Carrington Electrons Event Fluence'
ShieldingCurves['CarringtonElectronINTEGRALPowTabelated']['color'] = 'blue'

ShieldingCurves['ISS-AE9-mission']['label'] = '1 year dose ISS trapped Electrons'
ShieldingCurves['ISS-AP9-mission']['label'] = '1 year dose ISS trapped Protons'

ShieldingCurves['VAB-AE9-mission']['label'] = '1 month dose VAB trapped Electrons'
ShieldingCurves['VAB-AP9-mission']['label'] = '1 month dose VAB trapped Protons'


NumTiles = len(ShieldingCurves[Names[0]]['dose'])

## Plot shielding curves with error bars
x = np.linspace(0, 10, num=101, dtype=float, endpoint=True)
# print(x)


Ex = ShieldingCurves['Carrington-SEP-Expected-Int-With0']

for key, Curve in ShieldingCurves.items():
    # Fill the area between the expected and 2 sigma curves
    if "Plus" in key:
        plt.fill_between(x, Curve['dose'], Ex['dose'], color=Curve['color'], alpha=0.5, label=Curve['label'])
    elif "Minus" in key:
        plt.fill_between(x, Ex['dose'], Curve['dose'], color=Curve['color'], alpha=0.5, label=Curve['label'])
    else:
        plt.errorbar(x, Curve['dose'], Curve['error'], fmt='', markersize=2, capsize=2, label=Curve['label'], color=Curve['color'], linestyle='')
    

####### Plot 1kRad line #########
CriticalDose = [1 for i in x]
plt.plot(x, CriticalDose, color='k', linewidth=2, label='1 krad')
CriticalDose = [10 for i in x]
plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='10 krad')
# CriticalDose = [100 for i in x]
# plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='100 krad')
 """

plt.title("Ionising Dose Behind Shielding of Varying Thickness")
plt.xlabel("Aluminium Shielding Thickness [mm]")
plt.ylabel("Total Ionising Dose [krad]")
plt.yscale("log")
# plt.xscale("log")
# Increase the number of ticks on the y-axis
# plt.yticks(np.logspace(-4, 4, num=9))

plt.xlim(0.5, 10)
plt.ylim(1e-4, 2e+2)
plt.grid(which='both')

# Get handles and labels from the plot
handles, labels = plt.gca().get_legend_handles_labels()

# Define desired order of labels
desired_order = [
    '1 krad', 
    '10 krad', 
    'Carrington Proton + 2σ',
    'Carrington Proton Event Fluence',
    'Carrington Proton − 2σ',
    'Carrington Electrons Event Fluence',
    '1 year dose ISS trapped Electrons',
    '1 year dose ISS trapped Protons',
    '1 month dose VAB trapped Electrons',
    '1 month dose VAB trapped Protons'
]

# Reorder handles and labels according to the desired order
ordered_handles = [handles[labels.index(label)] for label in desired_order if label in labels]
ordered_labels = [label for label in desired_order if label in labels]

plt.legend(ordered_handles, ordered_labels)




plt.savefig(Path + "CompareGRASShieldose.pdf", format='pdf', bbox_inches="tight")
# plt.show()
