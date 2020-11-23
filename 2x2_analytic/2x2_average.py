# FYS3150 - Project 4, 2020.
# Simple program for the 2x2 lattice.
# Calculating average values and thermodynamic quantities
# using analytical expressions.

import numpy as np

J = 1            # Magnetic material constant.
k = 1            # Boltzmann constant.
T = 1            # Temperature.
B = 1/(k*T)      # Beta.

# Partition function Z.
Z = 4*np.cosh(8*B*J) + 12

# Mean energy <E>.
E = 32*J/Z * np.sinh(-8*B*J)

# Mean square energy <E^2>.
E2 = 256*J**2/Z * np.cosh(8*B*J)

# Energy variance sigma_E^2.
E_var = E2 - E**2

# Heat capacity Cv.
Cv = E_var/(k*T**2)

# Mean magnetization <M>.
M = 0

# Mean abs magnetization <|M|>.
M_abs = 8/Z * (np.exp(8*B*J) + 2)

# Mean square magnetization <M^2>.
M2 = 32/Z * (np.exp(8*B*J) + 1)

# Magnetization variance sigma_M^2.
M_var = M2 - M**2

# Magnetic susceptibility X.
X = M_var/(k*T)

# Alternative way, magnetic susceptibility:
M_var_alt = M2 - M_abs*M_abs
X_alt = M_var_alt/(k*T)

print("Z = {:.3g}".format(Z))
print("E = {:.3g}".format(E))
print("E2 = {:.3g}".format(E2))
print("E_var = {:.3g}".format(E_var))
print("Cv = {:.3g}".format(Cv))
print("M = {:.3g}".format(M))
print("M_abs = {:.3g}".format(M_abs))
print("M2 = {:.3g}".format(M2))
print("M_var = {:.3g}".format(M_var))
print("X = {:.3g}".format(X))
print("X_alt = {:.3g}".format(X_alt))
