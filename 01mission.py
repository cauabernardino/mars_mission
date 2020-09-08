import numpy as np
from toolbox import *


# For a circular orbit, r_dot = r . theta_dot (i_hat_theta)
R_DOT_NANO = R_LMO * THETA_DOT_LMO  # km/s


###################################################################

#### Task 1 ####

# Initial LMO orbit position
LMO_EA_T0 = np.array([OMEGA_LMO, I_LMO, THETA_LMO_T0])

LMO_EA_RATE = np.array([0, 0, THETA_DOT_LMO]) # Constant velocity, already in rad

LMO_R = np.array([R_LMO, 0, 0]).T  # Constant position vector


# Equivalent to h_hat = HN * n_hat
matrix_HN = EAtoDCM313(LMO_EA_T0)

# Mapping inverse - n_hat = NH * h_hat -> to find r_inertial 
matrix_NH = matrix_HN.T

r_inertial_T0 = np.matmul(matrix_NH, LMO_R)


# History of Euler Angles variation
EA_history = [LMO_EA_T0]

# History of positions in inertial frame
r_inertial = [r_inertial_T0]

# Initialization of integrations
x0 = LMO_EA_T0

for t in range(1, 451):
    x = x0 + t * EArate313(x0, LMO_EA_RATE)  
    EA_history.append(x)
    x0 = x

## Find r_inertial  in 450s

matrix_HN_450 = EAtoDCM313(EA_history[450])
matrix_NH_450 = matrix_HN_450.T
r_inertial_450 = np.matmul(matrix_NH_450, LMO_R.T)
v_inertial_450 = np.matmul(matrix_NH_450, LMO_EA_RATE.T)

print(v_inertial_450)
print(r_inertial_450)
