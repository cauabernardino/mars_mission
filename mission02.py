import numpy as np
from toolbox import *

"""
Task 3 - Sun-Pointing Reference Frame Orientation

Determine an analytic expressions for the sun pointing reference frame Rs by defining the DCM [RsN]

As the axis are fixed, there is no angular velocity.
"""

SUN_FRAME = np.array([OMEGA_SUN, I_SUN, 0])

matrix_RsN = np.round(EAtoDCM313(SUN_FRAME), 6)

print(f"Matrix RsN at 0s:\n{matrix_RsN}")

"""
Task 4 - Nadir-Pointing Reference Frame Orientation

Determine an analytic expression for the nadir pointing reference frame Rn by defining the DCM [RnN]

"""
## Matrix RnH
matrix_RnH = np.round(EAtoDCM313(NADIR_FRAME), 6)

## Matrix HN
# Initialization of LMO Euler Angles variation
EA_LMO_history = []

# LMO integration
x0 = LMO_EA_T0
step = 1
time = 330

for t in range(time):
    x = x0 + step * EArate313(x0, LMO_EA_RATE) 
    EA_LMO_history.append(x)
    x0 = x

## Find HN at 330s

matrix_HN = EAtoDCM313(EA_LMO_history[time - 1])

## Matrix RnN at 330s
matrix_RnN = np.round(np.matmul(matrix_RnH, matrix_HN), 6)

print(f"Matrix RnN at 330s:\n{matrix_RnN}")


## Angular velocity in inertial frame

matrix_NRn = matrix_RnN.T
v_inertial_nadir_330 = np.matmul(matrix_NRn, R_DOT_NANO.T)

print(f"Angular velocity in N frame at 330s:\n{v_inertial_nadir_330}")