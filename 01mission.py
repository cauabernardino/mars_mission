import numpy as np
from toolbox import *

"""
Task 1 - Orbit simulation

Find inertial position and velocity for LMO at 450s and GMO at 1150s

# For a circular orbit, R_DOT = R . THETA_DOT (i_hat_theta)
"""

#### LMO ####

# Initial LMO orbit position
LMO_EA_T0 = np.array([OMEGA_LMO, I_LMO, THETA_LMO_T0])
LMO_EA_RATE = np.array([0, 0, THETA_DOT_LMO]) # Constant velocity, in rad
LMO_R = np.array([R_LMO, 0, 0]).T  # Constant position vector
R_DOT_NANO = np.array([0, R_LMO * THETA_DOT_LMO, 0]) # Velocity in km/s

# Equivalent to h_hat = HN * n_hat
matrix_HN_LMO = EAtoDCM313(LMO_EA_T0)
# Mapping inverse - n_hat = NH * h_hat -> to find r_inertial 
matrix_NH_LMO = matrix_HN_LMO.T
r_inertial_LMO_T0 = np.matmul(matrix_NH_LMO, LMO_R)

# Initialization of LMO Euler Angles variation
EA_LMO_history = []

# LMO integration
x0 = LMO_EA_T0
step = 1

for t in range(450):
    x = x0 + step * EArate313(x0, LMO_EA_RATE) 
    EA_LMO_history.append(x)
    x0 = x

## Find LMO r and v at 450s

matrix_HN_LMO_450 = EAtoDCM313(EA_LMO_history[449])
matrix_NH_LMO_450 = matrix_HN_LMO_450.T
r_inertial_LMO_450 = np.matmul(matrix_NH_LMO_450, LMO_R.T)
v_inertial_LMO_450 = np.matmul(matrix_NH_LMO_450, R_DOT_NANO.T)

print(f"rLMO at 450s is: {r_inertial_LMO_450}")
print(f"vLMO at 450s is: {v_inertial_LMO_450}")


#### GMO ####

# Initial GMO orbit position
GMO_EA_T0 = np.array([OMEGA_GMO, I_GMO, THETA_GMO_T0])  # degrees
GMO_EA_RATE = np.array([0, 0, THETA_DOT_GMO])
GMO_R = np.array([R_GMO, 0, 0])
R_DOT_MOTHER = np.array([0, R_GMO * THETA_DOT_GMO, 0])  # Velocity in km/s

# Mapping through frames
matrix_HN_GMO = EAtoDCM313(GMO_EA_T0)
matrix_NH_GMO = matrix_HN_GMO.T
r_inertial_GMO_T0 = np.matmul(matrix_NH_GMO, GMO_R)

# Initialization of GMO Euler Angles variation
EA_GMO_history = []

# GMO integration
y0 = GMO_EA_T0
step = 1

for t in range(1150):
    y = y0 + step * np.rad2deg(GMO_EA_RATE)
    EA_GMO_history.append(y)
    y0 = y


## Find GMO r and v at 1150s

matrix_HN_GMO_1150 = EAtoDCM313(EA_GMO_history[1149])
matrix_NH_GMO_1150 = matrix_HN_GMO_1150.T
r_inertial_GMO_1150 = np.matmul(matrix_NH_GMO_1150, GMO_R.T)
v_inertial_GMO_1150 = np.matmul(matrix_NH_GMO_1150, R_DOT_MOTHER.T)

print(f"rGMO at 1150s is: {r_inertial_GMO_1150}")
print(f"vGMO at 1150s is: {v_inertial_GMO_1150}")


"""
Task 2 - Orbit Frame Orientation

Find the LMO DCM [HN] at 300s
"""

matrix_HN_LMO_300 = EAtoDCM313(EA_LMO_history[299])

print(f"The matrix HN at 300s is:\n{matrix_HN_LMO_300}")
