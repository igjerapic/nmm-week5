import os 
import glob
import re

import MDAnalysis as mda
import numpy as np
import beadspring as bsa
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from scipy.optimize import curve_fit

#def calculate_COM(atom_fragments):


def main():
    # Define the topology and trajectory files
    topology = 'c60_pair_with_bonds.data'

    traj_dir = "lmp_data"
    traj_files = glob.glob(os.path.join(traj_dir, "dump.data*"))
    traj_files = sorted(traj_files, key=lambda x : int(re.search(r'\d+', x).group()))

    u = bsa.setup_universe(topology, traj_files)

    N_FRAMES = u.trajectory.n_frames
    N_ATOMS_A = u.atoms.fragments[0].n_atoms
    N_ATOMS_B = u.atoms.fragments[-1].n_atoms
     
    # iniitiializations       
    c60a_COM = np.zeros((N_FRAMES, 3))
    c60b_COM = np.zeros((N_FRAMES, 3))

    c60a_forces = np.zeros((N_FRAMES, 3))
    c60b_forces = np.zeros((N_FRAMES, 3))

    # Loop over the frames and load the positions and forces
    for i, traj in enumerate(u.trajectory):    
        c60a  = u.atoms.fragments[0]
        c60b = u.atoms.fragments[-1]                      
        c60a_forces[i] = sum(c60a.forces.copy())
        c60b_forces[i] = sum(c60b.forces.copy())
        
        c60a_COM[i] = c60a.center_of_mass()
        c60b_COM[i] = c60b.center_of_mass()

    # vector pointing from c60a COM to c60b COM
    vect_connect = c60a_COM - c60b_COM

    # projection of forces on molecule a onto vector connecting of COMs
    distances = np.linalg.norm(vect_connect, axis=1)
    force_projection = np.einsum('ij,ij->i', c60a_forces, vect_connect ) / distances


    # Setting up bins
    min_dist = min(distances)
    max_dist = max(distances)
    num_bins = 50

    bin_edges = np.linspace(min_dist, max_dist, num_bins)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # averaging forces on molecule A
    avg_forces = np.zeros_like(bin_centers)
    for i in range(len(bin_centers)):
        mask = (distances >= bin_edges[i]) & (distances <= bin_edges[i + 1])
        avg_forces[i] = np.average(force_projection[mask])

    # removing nan values
    valid_indices = ~np.isnan(avg_forces)
    avg_forces_valid = avg_forces[valid_indices]
    bin_centers_valid= bin_centers[valid_indices]


    # determining Potential Mean Field 
    potential= np.zeros(len(bin_centers_valid))
    for i in range(len(bin_centers_valid) - 1, 0, -1):
        potential[i-1] = potential[i] + trapezoid(avg_forces_valid[i-1:i+1], 
                                             x= bin_centers_valid[i-1:i+1] )

    # fitting to LJ potential with offset
    def lj(r, epsilon, sigma):
        return 4 * epsilon * ( (sigma / r) ** 12 - (sigma / r) ** 6)

    bounds = ((0,0), (np.inf, np.inf))
    params, cov = curve_fit(lj, bin_centers_valid, potential,p0=(1,8), 
                            bounds = bounds)
    params_err = np.diag(cov)**0.5
    
    fit_vals = lj(bin_centers_valid, *params)


    # Plotting results
    plt.plot(bin_centers_valid, avg_forces_valid, "o:", color = "C1")
    plt.xlabel(r"Distance ($\AA$)")
    plt.ylabel("Force (eV / $\AA$)")
    plt.show()
    
    plt.plot(bin_centers_valid, potential, "o", label = "PMF")
    lj_label = f"LJ fit: $\epsilon$={params[0]:.2f}$\sigma={params[-1]:.2f} \pm$"
    plt.plot(bin_centers_valid, fit_vals, color = 'k' , )
    plt.xlabel(r"Distance ($\AA$)")
    plt.ylabel("Energy (eV)")
    plt.show()

if __name__ == '__main__':
    main()