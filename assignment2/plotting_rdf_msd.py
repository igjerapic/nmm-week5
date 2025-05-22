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
    epsilons = [0.5, 2.0, 5.0] 
    etas = [0.4, 0.8]

    for epsilon in epsilons:
        for eta in etas:
            file_name = f"epsilon{epsilon}_eta{eta}/" + data_file
            df = {}
            with open(file_name, 'rb') as f: 
                df = pkl.load(f)
            
            plt.figure(1)
            plt.plot(df["rdf_bincenters"], df["rdf"], label = f"$\eta=${eta}")

            plt.figure(2)
            plt.loglog(df["time_log"][1:], df["msd"], label=f"$\eta=${eta}")

        plt.figure(1)
        plt.xlabel(r"$r/\sigma$")
        plt.ylabel(r"$g(r)$")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"RDF_epsilon{epsilon}.svg", dpi=300)

        plt.figure(2)
        plt.xlabel(r"Time $\tau$")
        plt.ylabel(r"MSD")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"MSD_epsilon{epsilon}.svg", dpi=300)

        # Clearing Figures
        plt.figure(1)
        plt.clf()
        plt.figure(2)
        plt.clf()

    eta = 0.4
    epsilons = [0.5, 2.0, 5.0]
    for epsilon in epsilons:
        file_name = f"epsilon{epsilon}_eta{eta}/" + data_file
        df = {}
        with open(file_name, 'rb') as f: 
            df = pkl.load(f)
        
        plt.figure(1)
        plt.plot(df["rdf_bincenters"], df["rdf"], label = f"$\epsilon=${epsilon}")

        plt.figure(2)
        plt.loglog(df["time_log"][1:], df["msd"], label=f"$\epsilon=${epsilon}")

    plt.figure(1)
    plt.xlabel(r"$r/\sigma$")
    plt.ylabel(r"$g(r)$")
    plt.legend()
    plt.tight_layout()
    #plt.savefig("RDF_eta0.4.svg", dpi=300)

    plt.figure(2)
    plt.xlabel(r"Time $\tau$")
    plt.ylabel(r"MSD")
    plt.legend()
    plt.tight_layout()
    #plt.savefig("MSD_eta0.4.svg", dpi=300)
        
    plt.show()

if __name__=="__main__":
    main()