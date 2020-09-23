import numpy as np
from toolbox import *

"""
Task 1 - Orbit simulation

Find inertial position and velocity for LMO at 450s and GMO at 1150s

# For a circular orbit, R_DOT = R . THETA_DOT (i_hat_theta)
"""
#### LMO ####

# ## Find LMO r and v at 450s

matrix_HN_LMO_450 = orbit_integrator('LMO', 450, 1)
matrix_NH_LMO_450 = matrix_HN_LMO_450.T

r_inertial_LMO_450 = np.matmul(matrix_NH_LMO_450, LMO_R.T)
v_inertial_LMO_450 = np.matmul(matrix_NH_LMO_450, R_DOT_NANO.T)

print(f"rLMO at 450s is: {r_inertial_LMO_450}")
print(f"vLMO at 450s is: {v_inertial_LMO_450}")


#### GMO ####


## Find GMO r and v at 1150s

matrix_HN_GMO_1150 = orbit_integrator('GMO', 1150, 1)
matrix_NH_GMO_1150 = matrix_HN_GMO_1150.T

r_inertial_GMO_1150 = np.matmul(matrix_NH_GMO_1150, GMO_R.T)
v_inertial_GMO_1150 = np.matmul(matrix_NH_GMO_1150, R_DOT_MOTHER.T)

print(f"rGMO at 1150s is: {r_inertial_GMO_1150}")
print(f"vGMO at 1150s is: {v_inertial_GMO_1150}")


"""
Task 2 - Orbit Frame Orientation

Find the LMO DCM [HN] at 300s
"""

matrix_HN_LMO_300 = orbit_integrator('LMO', 300, 1)

print(f"The matrix HN at 300s is:\n{matrix_HN_LMO_300}")
