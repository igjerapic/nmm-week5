import numpy as np
import beadspring as bsa
import pickle as pkl

def main():
    # Define the topology and trajectory files

    topology = "final.data"
    traj_lin = "traj_lin.dat"
    traj_log = "traj_log.dat"

    u_lin = bsa.setup_universe(topology, traj_lin)
    u_log = bsa.setup_universe(topology, traj_log)

    N_FRAMES_LIN = u_lin.trajectory.n_frames
    N_FRAMES_LOG = u_log.trajectory.n_frames
    N_ATOMS = u_lin.atoms.n_atoms

    # Initialise the position and time arrays
    positions_lin = np.zeros((N_FRAMES_LIN, N_ATOMS, 3))
    time_lin = np.zeros(N_FRAMES_LIN)

    # Loop over the trajectory and load the positions
    for i,traj in enumerate(u_lin.trajectory):                          
        positions_lin[i] = u_lin.atoms.positions   
        time_lin[i] = u_lin.trajectory.ts.data['time']

    # Initialise the position and time arrays
    positions_log = np.zeros((N_FRAMES_LOG, N_ATOMS, 3))
    time_log = np.zeros(N_FRAMES_LOG)

    # Loop over the trajectory and load the positions
    for i,traj in enumerate(u_log.trajectory):                          
        positions_log[i] = u_log.atoms.positions   
        time_log[i] = u_log.trajectory.ts.data['time']
    
    # cubic box based on the dimensions of the simulation box
    box = bsa.setup_freud_box(u_lin.dimensions[0])

    # Computing RDF from final timesnap
    rdf_bincenters, rdf, r_min, r_peak = bsa.compute_rdf(positions_lin[-1], box, r_max=5.0, bins=50)

    # Computing MSD
    msd = bsa.compute_msd(positions_log)

    df_post_process = { "msd": msd,
                        "time_log": time_log,
                        "rdf_bincenters": rdf_bincenters,
                        "rdf": rdf,
                        }
    
    with open("post_process.pkl", "wb") as f:
        pkl.dump(df_post_process, f)


if __name__ == '__main__':
    main()
