import re, yaml


import pandas as pd
import matplotlib.pyplot as plt

from cycler import cycler
import lammps_logfile

plt.style.use('../../scripts/default.mplstyle')

plt.rcParams['axes.prop_cycle'] = plt.cycler(cycler(color = ['#332288', 
                                    '#CC6677',
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
    # reading non yaml LAMMPS log file
    log = lammps_logfile.File("log.lammps")
    df = pd.DataFrame(log.partial_logs[-1])
    df.to_pickle("thermo.pkl")
    
    keywords = [["Temp"], ["Press"], ['TotEng']]
    labels = ["Temp ($\epsilon/k_B)$", "Press ($\epsilon/\sigma^3$)", 'Total Energy ($\epsilon$)']
    for y, ylabel in zip(keywords, labels):
        fig = df.plot(x="Step", y=y, ylabel=ylabel)
        plt.tight_layout()
        plt.show()
if __name__ == '__main__':
    main()