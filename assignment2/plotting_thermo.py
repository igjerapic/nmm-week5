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
    
    eta = 0.4
    epsilons = [0.5, 2.0, 5.0]
    for epsilon in epsilons:
        file_name = f"epsilon{epsilon}_eta{eta}/" + data_file
        df = {}
        with open(file_name, 'rb') as f: 
            df = pkl.load(f)
        
        plt.figure(1)
        plt.plot(df["Time"], df["Temp"], label = f"$\\epsilon=${epsilon}")

        plt.figure(2)
        plt.plot(df["Time"], df["Press"], label=f"$\\epsilon=${epsilon}")

        plt.figure(3)
        plt.plot(df["Time"], df["TotEng"], label=f"$\\epsilon=${epsilon}")

    plt.figure(1)
    plt.xlabel(r"Time ($\tau$)")
    plt.ylabel(r"Temp $(\epsilon / k_B )$")
    plt.legend()
    plt.tight_layout()
    plt.savefig("temp_eta0.4.svg", dpi=300)
    
    plt.figure(2)
    plt.xlabel(r"Time ($\tau$)")
    plt.ylabel(r"Pressure ($\epsilon / \sigma^3$)")
    plt.legend(loc=(0.7, 0.6))
    plt.tight_layout()
    plt.savefig("press_eta0.4.svg", dpi=300)
   
    
    plt.figure(3)
    plt.xlabel(r"Time ($\tau$)")
    plt.ylabel(r"Total Energy ($\epsilon$)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("totEng_eta0.4.svg", dpi=300)
    #plt.show()

if __name__=="__main__":
    main()