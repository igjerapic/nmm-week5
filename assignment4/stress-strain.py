import re, yaml
import pickle as pkl

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import beadspring as bsa

INTEGRATION_TIME = 0.005

PLATEAU_STRAIN = 1.0 # from visual inspection, each stress-strain curve seems to plateau after staring of 1%

erates = np.array([0.001, 0.005, 0.01, 0.5, 0.1, 0.2, 0.4, 0.5, 0.8, 1.0])

#erates = np.array([1.0])
average_stress = np.zeros_like(erates)
viscosity = np.zeros_like(erates)

for i, erate in enumerate(erates):
    thermos_file = f"rate_{erate}/thermo.pkl" 
    df = pkl.load(open(thermos_file, "rb"))
    
    print(df.columns)
    stress = df["c_s4"].to_numpy()
    strain = df["Step"].to_numpy() * erate * INTEGRATION_TIME
    
    # plotting stress-strain curves and saving 
    plt.semilogx(strain, stress, ".")
    plt.title(f"Strain rate of {erate}")
    plt.show()
    plt.clf()

    # determining viscosity
    average_stress[i] = np.mean(stress[strain > PLATEAU_STRAIN])
viscosity = average_stress / erate

plt.loglog(erates, viscosity)
plt.show()
