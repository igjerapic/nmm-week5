import pickle as pkl
from cycler import cycler


import matplotlib.pyplot as plt
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
    with open("thermo.pkl", "rb") as f:
        df = pkl.load(f)
    
    keywords = [["Temp"], ["Press"], ['TotEng']]
    labels = ["Temp", "Press", 'Energy']

    for y, ylabel in zip(keywords, labels):
        fig = df.plot(x="Time", y=y, ylabel=ylabel, xlabel=r"Time $\tau$")
        plt.tight_layout()
        plt.show()
    # fig_temp = df.plot(x="Time", y="Temp", ylabel="Temp", figsize=(6,6))
    #plt.savefig('thermo_bondeng.png')
    # plt.show()

    # fig_Press = df.plot(x="Time", y="Press", ylabel="Press")
    # plt.show()

    # fig_energy = df.plot(x="Time", y=['KinEng', 'E_pair'], ylabel="Energy in Reduced Units")
    # plt.show()

if __name__ == '__main__':
    main()