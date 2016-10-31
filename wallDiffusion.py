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

# Consumption = - r_a * V
# V(r) = pi * L * (r^2 - r_0^2)
# r_a = k * C_a * d/dt (C_a)
# No idea if the above reaction is correct.


# Boundary conditions
# C_a = 0 @ r = R
# C_a = C_wall @ r = r_0

# This all boils down to :

# D_wall / r * d/dr (r * d/dr(C_a)) =  k * C_a * d/dt (C_a)

C_wall = 0 # The concentration at the wall. Defined from Jason's script.
C_final = 0 # The concentration at the nerve endings.
D_wall = 0 # Diffusion coefficient through the intestine wall
k = 0 # Rate coefficient of Tryp -> Ser.

# Now we just need to solve the differiential equation.
