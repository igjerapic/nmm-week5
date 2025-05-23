import re, yaml
import pandas as pd
import matplotlib.pyplot as plt

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

    df = pd.DataFrame(data=thermo[0]['data'], columns=thermo[0]['keywords'])
    
    df.to_pickle("thermo.pkl")
    
if __name__ == '__main__':
    main()