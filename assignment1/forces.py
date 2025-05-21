import os 
import glob

import MDAnalysis as mda
import numpy as np
import beadspring as bsa

#def calculate_COM(atom_fragments):


def main():
    # Define the topology and trajectory files
    topology = 'c60_pair_with_bonds.data'

    traj_dir = "lmp_data"
    traj_files = sorted(glob.glob(os.path.join(traj_dir, "dump.data*")))

    u = mda.Universe(topology, traj_files, format= "LAMMPSDUMP")

    N_FRAMES = u.trajectory.n_frames
    N_ATOMS_A = u.atoms.fragments[0].n_atoms
    N_ATOMS_B = u.atoms.fragments[-1].n_atoms
 
    # iniitiializations       
    c60a_positions = np.zeros((N_FRAMES, N_ATOMS_A, 3))
    c60b_positions = np.zeros((((N_FRAMES, N_ATOMS_B, 3))))

    c60a_forces = np.zeros_like(c60a_positions)
    c60b_forces = np.zeros_like(c60b_positions)

    # Loop over the frames and load the positions and forces
    for i,traj in enumerate(u.trajectory):    
        c60a  = u.atoms.fragments[0]
        c60b = u.atoms.fragments[-1]                      
        
        c60a_positions[i] = c60a.positions.copy()
        c60b_positions[i] = c60b.positions.copy()

        c60a_forces[i] = c60a.forces.copy()
        c60b_forces[i] = c60b.forces.copy()

    # total forces on molecules 
    c60a_force = np.sum(c60a_forces, axis=1)
    c60b_force = np.sum(c60b_forces, axis=1)

    # COM of molecules and connecting vector
    c60a_COM = np.sum(c60a_positions, axis=1) # allowed as all atoms have same mass
    c60b_COM = np.sum(c60b_positions, axis=1)

    vect_connect = c60a_COM - c60b_COM

    # projection of forces on molecule a onto connection of COMs
    distances = np.einsum('ij,ij->i', vect_connect, vect_connect)
    force_project = np.einsum('ij,ij->i', vect_connect, c60a_force ) / distances
    
    
    #TODO binning force on c60a over distances 


    #TODO fitting with either moors or LJ potential


    #TODO in.fullerene_CG

if __name__ == '__main__':
    main()