# If you are running this script on Habrok, don't forget to load Python2!
# First run: module load Python/2.7.18-GCCcore-12.2.0-bare

import random
from patch import patch

# Choose a relatively large lattice size keeping the runtime in mind
L = random.randint(20, 60)       # lattice size
L = 32
# Total number of particles
total_particles = L * L

# Split the total_particles into two species Nx and Ny
Nx = total_particles / 2 
Ny = total_particles - Nx

# Random patch counts per particle type
# x = random.randint(2, 4)       # x equidistant patches
# y = random.randint(2, 4)       # y equidistant patches
x = 2
y = 4



# Avoid large packing fractions, otherwise
# the system will be stuck in a crystalline phase
#vfrac = random.uniform(0, 0.74)
vfrac = 0.2

# Initialise patchy particle system
p = patch(vfrac=vfrac)
p.dim = 2
p.lattice = [L, L]              # square lattice of size L x L

# Build particles using 'ring' geometry

# ring: diam,N,1,2 = one large part with equatorial ring of N, types 1,2

# we fix the dimater to 2, set the radius of type 1 particles to 1.0  in OVITO
p.build(Nx, "ring", 2, x, 1, 2)
p.build(Ny, "ring", 2, y, 1, 2)

# Write output
p.write("init.data")
