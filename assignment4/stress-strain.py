import re, yaml
import pickle as pkl

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import beadspring as bsa

INTEGRATION_TIME = 0.005

PLATEAU_STRAIN = 1.0 # from visual inspection, each stress-strain curve seems to plateau after staring of 1%

erates = np.array([0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.4, 0.5, 0.8, 1.0])

#erates = np.array([1.0])
average_stress = np.zeros_like(erates)
average_stress_err = np.zeros_like(erates)
viscosity = np.zeros_like(erates)

for i, erate in enumerate(erates):
    thermos_file = f"rate_{erate}/thermo.pkl" 
    df = pkl.load(open(thermos_file, "rb"))
    
    print(df.columns)
    stress = df["c_s4"].to_numpy() / df["Volume"].to_numpy()
    strain = df["Step"].to_numpy() * erate * INTEGRATION_TIME
    
    # plotting stress-strain curves and saving 
    plt.semilogx(strain, stress, ".")
    

    # determining viscosity
    average_stress[i] = np.mean(stress[strain > PLATEAU_STRAIN])
    average_stress_err[i] = np.std(stress[strain > PLATEAU_STRAIN])
    plt.title(f"Strain rate of {erate}")
    plt.show()
    plt.clf()
viscosity = average_stress / erate
viscosity_err = average_stress_err / erate
plt.errorbar(erates, viscosity, viscosity_err)
plt.xscale("log")
plt.yscale("log")
#plt.loglog(erates, viscosity)
plt.show()
