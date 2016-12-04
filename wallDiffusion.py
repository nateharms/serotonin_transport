# Diffusion through the intenstine wall.

# I am going to model this as diffusion of substance through a porous material.

# Thickness of intestines goes from r_0 to R.

# Mass balance
# In - Out = Consumption.

# Fick's Law of Diffusion
# In = - D_wall * del_C @ r = r
# Out = - D_wall * del_C @ r = r + del_r

# Reaction
# tryptophan -> Serotonin

# conversion only occurs on the surface of the intestines
# Consumption = 0
# V(r) = pi * L * (r^2 - r_0^2)
# r_a = k * C_a * d/dt (C_a)
# No idea if the above reaction is correct.


# Boundary conditions
# dCdt = 0 @ r = r_0 + del_r
# C = C_wall @ r = r_0

# This all boils down to :

# D_wall / r * d/dr (r * d/dr(C_a)) =  0

# r d^2/dr^2 (C_a) + d/dr(C_a) = 0

# this solves out to c_1 * ln(r) + c_2

# B.C. are as fo


# The following will be defined from previous parts of the
C_wall = 0 # The concentration at the wall. Defined from Jason's script.
C_final = 0 # The concentration at the nerve endings.
D_wall = 0 # Diffusion coefficient through the intestine wall
k = 0 # Rate coefficient of Tryp -> Ser.

# Now we just need to solve the differiential equation.

import numpy as np
from scipy.integrate import ode
from numpy import log as ln

class WallDiffusion:
    def __init__(self, serWallConc, serDiffusion, serFinalConc):
        self.serWallConc = serWallConc
        self.serDiffusion = serDiffusion
        self.serFinalConcentration = serFinalconc

    def solveDiffusion(self):




def massBal(C, r):
    return c_1 * ln (r) + c_2

def solveDiffusion():
    # tbh I'm not 100% sure how to organize this part of the script. But I'll get it working soon enough
    func = massBal(C, r)

def pend(Y, t):
    C, C_prime = Y
    dCdt = [C, t]
