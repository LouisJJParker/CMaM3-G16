# Welcome to input parameters user interface
# To avoid errors, please input acceptable parameter values
# Save this file before run any linked further programs
# DEV: B-Beale, M-H-Fahrudin, S-Li, L-Parker, P-Wakley-Skinnarland

# String length: L = <POSITIVE INTEGER>
L = 1

# Total simulation time: final_time = <POSITIVE INTEGER>
final_time = 12000

# Time step: dt = <POSITIVE FLOAT>
dt = 0.1

# RK steps: RK_n = <POSITIVE INTEGER>
RK_n = 4

# Number of oscillators: n_b = <POSITIVE INTEGER>
n_b = 65

# Young’s modulus: k_young = <POSITIVE FLOAT>
k_young = 0.1

# Mass per unit length of the string: rho = <POSITIVE FLOAT>
rho = 400

# Non-linear coefficient: alpha = <POSITIVE FLOAT>
alpha = 0.25

# Amplitude of initial condition: pert_ini = <POSITIVE FLOAT>
pert_ini = 1.0

# Initial condition setting
    # "init_type = 1": single sine
    # "init_type = 2": half sine
    # "init_type = 3": parabola
init_type = 1
