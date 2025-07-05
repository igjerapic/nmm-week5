import os 
import glob
import re

import MDAnalysis as mda
import numpy as np
import beadspring as bsa
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from scipy.optimize import curve_fit
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

    c60a_positions = np.zeros((N_FRAMES, N_ATOMS_A, 3))
    rgsa = np.zeros(N_FRAMES)

    # Loop over the frames and load the positions and forces
    for i, traj in enumerate(u.trajectory):    
        c60a  = u.atoms.fragments[0]
        c60b = u.atoms.fragments[-1]

        _, eigvals = bsa.compute_gyration_tensor(c60a.positions.copy())
        rgsa[i] = bsa.calculate_rg2(*eigvals)
        c60a_positions[i] = c60a.positions.copy()

        c60a_forces[i] = sum(c60a.forces.copy())
        c60b_forces[i] = sum(c60b.forces.copy())
        c60a_COM[i] = c60a.center_of_mass()
        c60b_COM[i] = c60b.center_of_mass()

    # vector pointing from c60a COM to c60b COM
    vect_connect = c60a_COM - c60b_COM

    # projection of forces on molecule a onto vector connecting of COMs
    distances = np.linalg.norm(vect_connect, axis=1)
    force_projection = np.einsum('ij,ij->i', c60a_forces, vect_connect ) / distances


    # radius of gyration of fullerne atoms
    print(np.mean(rgsa)**0.5 * 2)
    print(rgsa[0]**0.5 * 2)

    print("COMa", c60a_COM[0])
    print("COMb", c60b_COM[0])


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

    # fitting to LJ potential
    def lj(r, epsilon, sigma):
        return 4 * epsilon * ( (sigma / r) ** 12 - (sigma / r) ** 6)
        
    def morse(r, D_e, k, r_0):
        a = (k / (2 *  D_e))**0.5
        return D_e * (np.exp( - 2 * a * (r - r_0))- 2 * np.exp(-a * (r - r_0)))

    bounds = ((0,0), (np.inf, np.inf))
    params_lj, cov_lj = curve_fit(lj, bin_centers_valid, potential,p0=(1,8), 
                            bounds = bounds)
    params_lj_err = np.diag(cov_lj)**0.5
    
    fit_vals_lj = lj(bin_centers_valid, *params_lj)

    bounds = ((0,0,0), (np.inf, np.inf, np.inf))
    params_morse, cov_morse = curve_fit(morse, bin_centers_valid, potential,p0=(0.5,1.0, 10.), 
                            bounds = bounds)
    params_morse_err = np.diag(cov_morse)**0.5
    
    fit_vals_morse = morse(bin_centers_valid, 0.5, 1.0, 10)
    fit_vals_morse = morse(bin_centers_valid, *params_morse)
    # Plotting results
    plt.plot(bin_centers_valid, avg_forces_valid, "o:", color = "C1")
    plt.xlabel(r"Distance ($\AA$)")
    plt.ylabel("Force (eV / $\AA$)")
    plt.tight_layout()
    # plt.savefig("figs/averged_force.png", dpi = 300)
    # plt.show()
    plt.clf()


    plt.plot(bin_centers_valid, potential, "o", label = "PMF")
    lj_label = f"LJ fit:\n\t$\epsilon={params_lj[0]:.2f} \pm {params_lj_err[0]:.2f}$\n\t$\sigma={params_lj[-1]:.2f} \pm {params_lj_err[-1]:.2f}$"
    plt.plot(bin_centers_valid, fit_vals_lj, color = 'k' , linewidth = 3, label = lj_label)

    morse_label = f"Morse fit:\n\t$D_e = {params_morse[0]:.2f} \pm {params_morse_err[0]:.2f}$\n\t" \
                + f"$k={params_morse[1]:.2f} \pm {params_morse_err[1]:.2f}$\n\t" \
                + f"$r_0={params_morse[-1]:.2f} \pm {params_morse_err[-1]:.2f}$"
    plt.plot(bin_centers_valid, fit_vals_morse, color = 'grey' , linestyle = "--", linewidth = 3, label = morse_label)
    plt.legend()
    plt.xlabel(r"Distance ($\AA$)")
    plt.ylabel("Energy (eV)")
    plt.tight_layout()
    # plt.savefig("figs/PMF.png", dpi = 300)
    # plt.show()

if __name__ == '__main__':
    main()