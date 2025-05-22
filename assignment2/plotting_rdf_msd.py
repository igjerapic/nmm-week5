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

    data_file = "post_process.pkl"
    # determining for 
    eta = 0.45
    temps = [0.5, 1.0, 2.0]

    for T in temps:
        file_name = f"temp{T}_eta{eta}/" + data_file
        df = {}
        with open(file_name, 'rb') as f: 
            df = pkl.load(f)
        
        plt.figure(1)
        plt.plot(df["rdf_bincenters"], df["rdf"], label = f"T{T}")

        plt.figure(2)
        plt.loglog(df["time_log"][1:], df["msd"], label=f"T{T}")

    plt.figure(1)
    plt.xlabel(r"$r/\sigma$")
    plt.ylabel(r"$g(r)$")
    plt.legend()
    plt.tight_layout()
    plt.savefig("RDF_eta0.45.svg", dpi=300)

    plt.figure(2)
    plt.xlabel(r"Time $\tau$")
    plt.ylabel(r"MSD")
    plt.legend()
    plt.tight_layout()
    plt.savefig("MSD_eta0.45.svg", dpi=300)

    # Clearing Figures
    plt.figure(1)
    plt.clf()
    plt.figure(2)
    plt.clf()

    temp = 1.0
    etas = [0.5, 0.6, 0.7]
    for eta in etas:
        file_name = f"temp{temp}_eta{eta}/" + data_file
        df = {}
        with open(file_name, 'rb') as f: 
            df = pkl.load(f)
        
        plt.figure(1)
        plt.plot(df["rdf_bincenters"], df["rdf"], label = f"eta{eta}")

        plt.figure(2)
        plt.loglog(df["time_log"][1:], df["msd"], label=f"eta{eta}")

    plt.figure(1)
    plt.xlabel(r"$r/\sigma$")
    plt.ylabel(r"$g(r)$")
    plt.legend()
    plt.tight_layout()
    plt.savefig("RDF_temp1.0.svg", dpi=300)

    plt.figure(2)
    plt.xlabel(r"Time $\tau$")
    plt.ylabel(r"MSD")
    plt.legend()
    plt.tight_layout()
    plt.savefig("MSD_temp1.0.svg", dpi=300)
    
   # plt.show()

if __name__=="__main__":
    main()