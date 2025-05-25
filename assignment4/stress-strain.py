import re, yaml
import pickle as pkl
from cycler import cycler

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


plt.style.use('../scripts/default.mplstyle')

plt.rcParams['axes.prop_cycle'] = plt.cycler(cycler(color = ['#CC6677', 
                                    '#332288', 
                                    '#88CCEE',
                                    '#DDCC77', 
                                    '#117733', 
                                    '#882255', 
                                    '#44AA99', 
                                    '#999933', 
                                    '#AA4499',
                                    '#DDDDDD'
                                ]))

INTEGRATION_TIME = 0.005
PLATEAU_STRAIN = 1.0 # from visual inspection, each stress-strain curve seems to plateau after staring of 1%

erates = np.array([0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.4, 0.5, 0.8, 1.0])

#erates = np.array([1.0])
average_stress = np.zeros_like(erates)
average_stress_err = np.zeros_like(erates)

for i, erate in enumerate(erates):
    thermos_file = f"rate_{erate}/thermo.pkl" 
    df = pkl.load(open(thermos_file, "rb"))
    
    # determining stress and strain; LAMMPS saves stress in units of Pressure * volume
    stress = df["c_s4"].to_numpy() / df["Volume"].to_numpy()
    strain = df["Step"].to_numpy() * erate * INTEGRATION_TIME
    
    # determining viscosity
    average_stress[i] = np.mean(stress[strain > PLATEAU_STRAIN])
    average_stress_err[i] = np.std(stress[strain > PLATEAU_STRAIN])

    # plotting stress-strain curves and saving 
    plt.semilogx(strain[1:], stress[1:], ".:", label = f"erate={erate}")
    
# plotting stess-strain curves
plt.xlabel("Strain (%)")
plt.ylabel(r"Stress $(\epsilon /\sigma^3)$")
plt.legend(ncol = 2, fontsize = 9,columnspacing=0.8)
plt.tight_layout()
#plt.show()
plt.clf()


# Plotting Flow curve and Fitting power law
viscosity = average_stress / erates
viscosity_err = average_stress_err / erates


# power law for erate
power_law = lambda x, m , a : a * x** m 
params, cov = curve_fit(power_law, erates, viscosity)
params_err = np.diagonal(cov) ** 0.5

fit_vals = power_law(erates, *params)

print("No sigma:")
print(params)
print(params_err)

# weighted average for viscosity for lower bound of zero shear viscosity

print(np.mean(viscosity))
# plotting 
plt.errorbar(erates, viscosity, viscosity_err, linestyle = ":", marker=".", capsize=2)
plt.plot(erates, fit_vals, "grey", linestyle = ":")

plt.xlabel(r"Strain rate (1/$\tau$)")
plt.ylabel(r"Viscosity ($\epsilon\tau/\sigma^3$)")
plt.xscale("log")
plt.yscale("log")
plt.tight_layout()
plt.show()
