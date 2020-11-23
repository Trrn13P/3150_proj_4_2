# FYS3150 - Project 4, 2020.
# Ising Model for LxL lattice.
# Using periodic boundary conditions.

import numpy as np
import random

L = 4            # Number of spins in each dimension.
J = 1            # Magnetic material constant.
k = 1            # Boltzmann constant.
T = 1            # Temperature.
N = L*L          # Number of spins in LxL lattice.
beta = 1/(k*T)   # Dimension constant for Boltzmann factors.
MCS  = int(1e5)  # Monte Carlo Samples.

# Initial microstate: All spins pointing up, ground state:
lattice_0 = np.ones((L,L))   # LxL lattice with N = LxL spins.

# You can also initiate with random spins:
# lattice_0 = [[-1,1,1,1],[1,-1,-1,1],[-1,1,1,1],[-1,1,-1,-1]]

def PBC(k,limit,add):
    # Periodic boundary conditions.
    return (k + limit + add) % limit     # Indexing.

def energy(lattice):
    # Calculating energy for a given microstate in the lattice.
    energy = 0
    for i in range(L):
        for j in range(L):
            Eij = -J*(lattice[i][j]*(lattice[PBC(i,L,-1)][j] + lattice[i][PBC(j,L,-1)]))
            energy = energy + Eij
    return energy

def metropolis(MCS, lattice):
    # Metropolis algorithm, using Monte Carlo sampling.
    E = energy(lattice)
    for i in range(MCS):      # Monte Carlo cycles.
        for j in range(L):
            for k in range(L):
                r1 = random.uniform(0,1)         # Random float between 0 and 1.
                r2 = random.randint(0,L-1)       # Random integer between 0 and L-1.
                r3 = random.randint(0,L-1)       # Random integer between 0 and L-1.
                Ei = energy(lattice)             # Initial energy.
                lattice[r2][r3] = lattice[r2][r3]*(-1)   # Flip spin, random index.
                Ej = energy(lattice)             # New energy.
                dE = Ej - Ei                     # Energy difference.
                if dE <= 0:
                    E = E + dE            # Accept new configuration.
                else:
                    w = np.exp(-beta*dE)  # Relative probability, Boltzmann factor.
                    if r1 <= w:
                        E = E + dE        # Accept new configuration.
                    else:
                        # Deny new configuration, flip spin back.
                        lattice[r2][r3] = lattice[r2][r3]*(-1)
                        E = E + Ei
    return E, lattice

E, lattice = metropolis(MCS, lattice_0)
E_mean = E/MCS/N

print("Mean energy of the system, <E> = {:.4g}".format(E_mean))
print("Final microstate configuration:")
print(lattice)
