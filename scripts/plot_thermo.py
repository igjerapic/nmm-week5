import pickle as pkl

#import pandas as pd
import matplotlib.pyplot as plt

def main():
    with open("thermo.pkl", "rb") as f:
        df = pkl.load(f)
    
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