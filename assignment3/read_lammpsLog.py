import re, yaml
import pandas as pd
import matplotlib.pyplot as plt

from cycler import cycler
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
    try:
        from yaml import CSafeLoader as Loader
    except ImportError:
        from yaml import SafeLoader as Loader

    docs = ""
    with open("log.lammps") as f:
        for line in f:
            m = re.search(r"^(keywords:.*$|data:$|---$|\.\.\.$|  - \[.*\]$)", line)
            if m: docs += m.group(0) + '\n'

    thermo = list(yaml.load_all(docs, Loader=Loader))
    print(thermo)
    df = pd.DataFrame(data=thermo['data'], columns=thermo['keywords'])
    
    df.to_pickle("thermo.pkl")
    keywords = [["Temp"], ["Press"], ['KinEng', 'E_pair']]
    labels = ["Temp", "Press", 'Energy']

    for y, ylabel in zip(keywords, labels):
        fig = df.plot(x="Time", y=y, ylabel=ylabel)
        plt.tight_layout()
        plt.show()
    fig_temp = df.plot(x="Time", y="Temp", ylabel="Temp", figsize=(6,6))
    #plt.savefig('thermo_bondeng.png')
    plt.show()

    fig_Press = df.plot(x="Time", y="Press", ylabel="Press")
    plt.show()

    fig_energy = df.plot(x="Time", y=['KinEng', 'E_pair'], ylabel="Energy in Reduced Units")
    plt.show()
if __name__ == '__main__':
    main()