import os
import pickle as pkl

import numpy as np
import freud 

import beadspring as bsa
import MDAnalysis as mda
import matplotlib.pyplot as plt

def find_clusters(positions, box, CUTOFF=1.5):
    wrapped_coords = box.wrap(positions)
    system = freud.locality.AABBQuery(box,  wrapped_coords)


    cl = freud.cluster.Cluster()
    clp = freud.cluster.ClusterProperties()
    cl.compute(system, neighbors={"r_max": CUTOFF})
    clp.compute(system, cl.cluster_idx)
    clusters = {}
    for i in range(cl.num_clusters):
        clusters[i] = system.points[cl.cluster_keys[i]]


    return cl, clp, clusters

def plot_clusters(cl, system):
    fig, ax = plt.subplots(1, 1, figsize=(9, 6))
    for cluster_id in range(cl.num_clusters):
        cluster_system = freud.AABBQuery(
            system.box, system.points[cl.cluster_keys[cluster_id]]
        )
        cluster_system.plot(ax=ax, s=10, label=f"Cluster {cluster_id}")
        print(
            f"There are {len(cl.cluster_keys[cluster_id])} points in cluster {cluster_id}."
        )
        plt.show()

def process(dir):
    topology = dir + "/final.data"
    trajectory  = dir + "/traj.dat"

    u = bsa.setup_universe(topology, trajectory)

    box = bsa.setup_freud_box(u.dimensions[0], dimensions=2)

    # only interested in final frame 
    u.trajectory[-1]

    positions = u.atoms.positions.copy()

    cl, clp, clusters = find_clusters(positions, box, CUTOFF=1.5)
    
    num_clusters = cl.num_clusters
    rgs = clp.radii_of_gyration
    gyrations = clp.gyrations
    eigvals = np.linalg.eigvals(gyrations)
    sizes = clp.sizes

    rhs = np.zeros_like(rgs)
    k2s = np.zeros_like(rgs)
    cs = np.zeros_like(rgs)
    for i in range(num_clusters):
        rhs[i] = bsa.calculate_hydrodynamic_radius(clusters[i])
        k2s[i] = bsa.calculate_shape_anisotropy(*np.sort(eigvals[i]))
        cs[i] = bsa.calculate_acylindricity(*np.sort(eigvals[i])[1:])

    data = {'num_clusters' : num_clusters,
            'rgs' : rgs,
            'sizes': sizes,
            'clusters': clusters,
            'rhs': rhs,
            'cs': cs
            }
    
    
    mask = np.zeros_like(rgs, dtype=int)
    weights = sizes 


    for i in range(num_clusters):
        if len(clusters[i]) > 8:
            mask[i] = 1
    
    print(np.mean(cs[mask]))
    print(np.average(cs, weights = sizes))
    print(np.mean(rgs)/ np.mean(rhs))
    return data

def main():
    Nxs = [0.25, 0.5]
    Ts = [0.5, 2.0]

    for Nx in Nxs:
        for T in Ts:
            dir_name = f"Nx{Nx}_T{T}"
            print(dir_name)
            data = process(dir_name)
            with open(dir_name +'/clusters.pkl', 'wb') as f:
                pkl.dump(data, f)


if __name__ == "__main__":
    main()
