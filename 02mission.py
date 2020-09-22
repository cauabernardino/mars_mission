import numpy as np
from toolbox import *

"""
Task 3 - Sun-Pointing Reference Frame Orientation

Determine an analytic expressions for the sun pointing reference frame Rs by defining the DCM [RsN]

As the axis are fixed, there is no angular velocity.
"""

SUN_FRAME = np.array([OMEGA_SUN, I_SUN, 0])

matrix_RsN = EAtoDCM313(SUN_FRAME)

print(matrix_RsN)

"""
Task 4 - Nadir-Pointing Reference Frame Orientation

Determine an analytic expression for the nadir pointing reference frame Rn by defining the DCM [RnN]

It should be RnH * HN
"""