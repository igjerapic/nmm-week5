# nmm-week5

In celebration of spherical cows, this week introduces Coarse-Grained MD, where the beads of the simulation do not necessarily correspond to individual atoms. These models lose predictive power but gain tremendous flexibility. In this week's assignment, you will focus on separate particles with complex interatomic interactions, able to represent any colloidal material from Quantum Dots to Yogurt!

## Assignment 1 - What is a particle?

The key for building new CG models is to calibrate the definition and interactions of your CG beads, either with bottom-up (e.g. from AA simulations) or top-down (e.g. from experimental results) approaches. For simple spherical particles, a common approach is to calculate the [Potential of Mean Force (PMF)](https://en.wikipedia.org/wiki/Potential_of_mean_force) between two "particles". Let's do that for fullerene in vacuum!

### Instructions

1a. After visualizing the c60-bonds.data file, run the `in.fullerene` simulation. Create an empty directory `lmp_data/` in the same directory as the script.

(i) What is the role of the looping structure in this simulation?

(ii) From the thermodynamic information in the log file/the trajectory, extract two plots: the force and the potential energy between the two molecules as a function of the distance between their centers of mass.

> **Hints for the PMF calculation:** We save the forces acting on every atom (fx, fy, fz). Identify different molecules in the simulation using MDAnalysis (what does `u.atoms.fragments` do?). You can access the forces on particles just like positions: `u.atoms.forces`. Sum the forces over all atoms in each molecule to get the total force vector on molecule 1 and molecule 2. Compute the centre of mass of each molecule, then the vector connecting them. Project the force on one molecule onto this vector (why?). Store the distance and projected force for each frame. Then average the forces in distance bins and integrate using `scipy.integrate` (we let you decide which function to use) to get the PMF.

1b. Time to coarse-grain!

(i) Fit your potential energy plot with an interatomic potential function of your choice and report the values of the fitting parameters. 

(ii) How would you now define a simulation of coarse-grained particles interacting with your effective potential? Specify all the information you would need to set it up.

## Assignment 2 - Just Jamming

Not bound to the periodic table anymore! We can now design "atoms" with arbitrary interactions, begging the question: what structures will they form based on their interaction potentials? In this assignment, let's discover the phases of the well-known [Lennard-Jones fluid](https://en.wikipedia.org/wiki/Lennard-Jones_potential).

### Instructions

2a. Take the new `in.3dlj` script, using a cutoff of 2.5 for the Lennard-Jones potential. We placed question marks for the parameters you need to change. If you try to run the script as is, it will throw an error. Look for comments with `!!!!!!!!!!` to locate those parameters and adjust accordingly.

(i) Fix the temperature to 1.0, and vary the particle density and interaction strength. Report snapshots and quantitative metrics (on structure and/or dynamics) for different phases of the system, including at least a liquid phase and a crystalline phase.

> **Pro tip:** Equilibrate long enough, and stay far away from the phase boundaries unless you want to wait for a long time.

## Assignment 3 - Patch up

Not all colloidal particles have spherical, isotropic interactions! Directional interactions are often present as a result of internal chemical structure or functionalization. In this assignment, learn how directional interactions can be implemented by adding attractive patches on repulsive cores.

### Instructions 

3a. Explore the phase diagrams of patchy particle systems in 2D.

(i) Build at least two initial data files with different compositions of patchy particles, each particle having 2, 3, or 4 patches. To do so read the `create_patch.py` and replace the random numbers in the script with your system parameters. More comments in the script. The script uses a Python tool created for LAMMPS called [pizza.py](https://lammps.github.io/pizza/). It is written in Python2, so make sure to load it from the modules stack on Habrok. Also make sure to run `create_patch.py` in the same directory as the two module files provided: `patch.py` and `data.py` 

(ii) For each composition, use the `in.patchy` file and vary thermodynamics conditions (temperature, packing fraction) to report at least two different phases. Justify your choice of composition, equilibration time, and quantitative metrics to define the structures. 

> **Pro tip:** You can search the literature to guide your choices. Or just run random simulations and see what sticks, who are we to judge how you want to spend your time.

## Assignment 4 - Look Ma, no hands!

Equilibrium phase diagrams are fun, but colloids are very interesting non-Newtonian fluids! Let's check their rheology in this assignment. 

### Instructions

4a. Start from an equilibrated structure obtained with the `in.3dlj` script (see assignment 2) with 2048 particles, temperature 1.0, packing fraction 0.58, $\varepsilon = 1.0$.


(i) Run the simulation `in.shear`, and obtain the viscosity from the resulting stress curve. The simulations runs in 15 minutes at shear rate 1.0 with the above parameters. What is the algorithm used to stabilize the temperature in the system?

> **Hints for viscosity calculation:** Thermodynamic output is stored in YAML format, which can be conveniently loaded into a pandas dataframe. See the [documentation](https://docs.lammps.org/Howto_structured_data.html#yaml-format-thermo-style-or-dump-style-output) for more details on YAML-formatted thermo output. You can use the code snippet below to parse the `log.lammps` file. Once loaded, search for the column labelled `c_s4`, check the corresponding `compute` command in the LAMMPS input script to explain what this quantity represents. Plot the stress-strain curves (how can you obtain the values of stress at a given strain value and simulation time step?). The stress typically rises sharply at first (elastic response), reaches a maximum (the overshoot), and then levels off (steady flow). Extract the average value of the stress in this plateau region. Then, divide this average stress by the applied shear rate. This gives you the shear viscosity
$$\eta = \frac{\langle \sigma \rangle}{\dot{\gamma}}$$ where $\langle \sigma \rangle$ is the mean steady-state stress and $\dot{\gamma}$ is the shear rate.

```python
import re, yaml
import pandas as pd
try:
    from yaml import CSafeLoader as Loader
except ImportError:
    from yaml import SafeLoader as Loader

docs = ""
with open('log.lammps') as f:
    for line in f:
        m = re.search(r"^(keywords:.*$|data:$|---$|\.\.\.$|  - \[.*\]$)", line)
        if m:
            docs += m.group(0) + '\n'

thermo = list(yaml.load_all(docs, Loader=Loader))
df = pd.DataFrame(data=thermo[1]['data'], columns=thermo[1]['keywords'])
```




(ii) Repeat the simulation for smaller shear rates. Following the same procedure as before, make a plot of viscosity as a function of shear rate (in log-log scale). Can you fit a power law to it? Can you determine a lower limit for the zero shear viscosity? 

> **Pro tip:** Change shear rate logarithmically. Be mindful of the resources needed: anything lower than a shear rate of 0.001 will take at least several hours to run.

**4b. (OPTIONAL, HARD)** Calculate the zero shear viscosity. This can be done in a couple of ways [(or more, none easy)](https://docs.lammps.org/Howto_viscosity.html):

(i) Keep doing what you were doing in 4a, but MUUUCH slower. Arm yourself with patience, and good luck!

(ii) Use the Green-Kubo formalism to extract the viscosity from the ensemble average of the auto-correlation of the stress tensor. Good luck!


