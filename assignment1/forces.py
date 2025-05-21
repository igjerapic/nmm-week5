import os 
import glob
import re

import MDAnalysis as mda
import numpy as np
import beadspring as bsa
import matplotlib.pyplot as plt

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
    force_project = np.einsum('ij,ij->i', vect_connect, c60b_forces ) / distances


    stepsizes = [distances[i] - distances[i - 1] for i in range(1, len(distances))]
    steps = np.arange(1, len(distances) + 1)

    #TODO binning distances and then averageing force on c60a over distances 
    bin_count_large = 20
    bin_count_small = 10
    center = len(distances)//2

    bin_edges_large = np.linspace(min(distances[:center]), max(distances[:center]),
                                  bin_count_large)
    bin_centers_large = 0.5 * (bin_edges_large[:-1] + bin_edges_large[1:])


    bin_edges_small = np.linspace(min(distances[center - 5:]), max(distances[center - 5:]), 
                                  bin_count_small)
    bin_centers_small = 0.5 * (bin_edges_small[:-1] + bin_edges_small[1:])

    bin_centers = np.concatenate((bin_centers_large, bin_centers_small))
    bin_edges = np.concatenate((bin_edges_large, bin_edges_small))

    # averaging
    avg_forces = np.zeros_like(bin_centers)
    for i in range(len(bin_centers)):
        mask = (distances >= bin_edges[i]) & (distances <= bin_edges[i + 1])
        avg_forces[i] = np.average(force_project[mask])
    #TODO fitting with either moors or LJ potential


    #TODO in.fullerene_CG

if __name__ == '__main__':
    main()