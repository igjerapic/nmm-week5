import pickle as pkl
import os
from cycler import cycler
from collections import Counter


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



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

def main():
    # change working director to that of file
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    
    Nxs = [0.25, 0.5]
    temps = [0.5, 2.0]
    avg_sizes = np.zeros((2,2)) # average sizes with shape (Nxs, temps)
    for i, Nx in enumerate(Nxs):
        cluster_sizes = {}
        for j, T in enumerate(temps):
            file= f"Nx{Nx}_T{T}_clusters.txt"
            df = np.loadtxt(file, comments="#").T

            cluster_sizes[f"{T}"] = df[1] # cluster sizes are stored in second column
            N_ATOMS = sum(cluster_sizes[f"{T}"])

            # determining weighted avg, where weight is size of cluster
            avg_sizes[i,j] = np.sum(np.power(cluster_sizes[f"{T}"], 2)) / N_ATOMS

    # plotting
    barWidth = 0.2
    fig = plt.subplots() 
    # Subbars for avg of temp 0.5
    br1 = np.arange(len(avg_sizes[:, 0])) 
    br2 = [x + barWidth for x in br1]

    plt.bar(br1, avg_sizes[:,0], width = barWidth, label ='$T=$0.5') 
    plt.bar(br2, avg_sizes[:, 1], width=barWidth, label ='$T=$2.0')
    plt.ylabel('Weighted Average Cluster Size') 
    plt.xticks([r + 0.5 * barWidth for r in range(len(avg_sizes[0]))], 
            ['Nx0.25', 'Nx0.5'])
    plt.tight_layout()
    plt.legend()
    plt.savefig("avg_clusters.svg")
    #plt.show() 

if __name__=="__main__":
    main()