import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import os
from cycler import cycler


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

    data_file = "thermo.pkl"
    
    Nxs = [0.25, 0.5]
    temps = [0.5, 2.0]
    for Nx in Nxs:
        for T in temps:
            file_name = f"Nx{Nx}_T{T}/" + data_file
            df = {}
            with open(file_name, 'rb') as f: 
                df = pkl.load(f)
            
            plt.figure(1)
            plt.plot(df["Step"], df["Temp"], label = f"$T=${T}")

            plt.figure(2)
            plt.plot(df["Step"], df["Press"], label=f"$T=${T}")

            plt.figure(3)
            plt.plot(df["Step"], df["TotEng"], label=f"$T=${T}")

        file_label = f"Nx{Nx}"
        plt.figure(1)
        plt.xlabel(r"Step")
        plt.ylabel(r"Temp $(\epsilon / k_B )$")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"temp_{file_label}.svg", dpi=300)
        plt.clf()

        plt.figure(2)
        plt.xlabel(r"Step")
        plt.ylabel(r"Pressure ($\epsilon / \sigma^3$)")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"press_{file_label}.svg", dpi=300)
        plt.clf()
        
        plt.figure(3)
        plt.xlabel(r"Step")
        plt.ylabel(r"Total Energy ($\epsilon$)")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"totEng_{file_label}.svg", dpi=300)
        plt.clf()
        #plt.show()

if __name__=="__main__":
    main()