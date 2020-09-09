import numpy as np
from numpy import sin, cos, arctan, arccos

"""
Mission scenario and constants

"""

### Mars
# Mars rotational period is 1 day and 37 minutes
R_MARS = 3396.19 # km - Radius
MI_MARS = 42828.3 # km^3/s^2 - Gravity constant

### Nano-sat (LMO) - Circular
h = 400 # Altitude
THETA_DOT_LMO = 0.000884797 # rad/sec - Constant orbit rate
R_LMO = R_MARS + h # km - Orbit radius

### Mothercraft (GMO) - Circular
# Has the sabe rotational period than Mars
R_GMO = 20424.2 # km - Orbit radius
THETA_DOT_GMO = 0.0000709003 # rad/sec - Constant orbit rate

##### Initial conditions #####

### Orbits
OMEGA_LMO = 20
I_LMO = 30
THETA_LMO_T0 = 60

OMEGA_GMO = 0
I_GMO = 0
THETA_GMO_T0 = 250

### Nano-sat
sigma_BN_T0 = np.array([0.3, -0.4, 0.5]).T  # Initial MRP
omega_BN_T0 = np.deg2rad([1.00, 1.75, -2.20]) # rad/s - Initial body angular velocity

# Intertia Tensor
I_nano = np.array([10, 5, 7.5]) * np.eye(3)  # kg . m^2


"""
AUXILIARY FUNCTIONS

"""

def EAtoDCM313(t, radian=False):
    """Insert 3-1-3 Euler angles in degrees and returns equivalent DCM"""
    if not radian:
        t = np.deg2rad(t)
    
    a11 = cos(t[2])*cos(t[0]) - sin(t[2])*cos(t[1])*sin(t[0])
    a12 = cos(t[2])*sin(t[0]) + sin(t[2])*cos(t[1])*cos(t[0])
    a13 = sin(t[2])*sin(t[1])
    a21 = -sin(t[2])*cos(t[0]) - cos(t[2])*cos(t[1])*sin(t[0])
    a22 = -sin(t[2])*sin(t[0]) + cos(t[2])*cos(t[1])*cos(t[0])
    a23 = cos(t[2])*sin(t[1])
    a31 = sin(t[1])*sin(t[0])
    a32 = -sin(t[1])*cos(t[0])
    a33 = cos(t[1])

    DCM = np.array([[a11, a12, a13],
                    [a21, a22, a23],
                    [a31, a32, a33]])

    return DCM 


def DCMtoEA313(DCM, radian=False):
    """
    Insert DCM and returns Euler Angles.

    radian: True, returns radians. False, returns degrees
    """
    t1 = arctan(DCM[2][0]/(-DCM[2][1]))
    t2 = arccos(DCM[2][2])
    t3 = arctan(DCM[0][2]/DCM[1][2])

    EA = np.array([t1, t2, t3])

    if radian:
        return EA
    else:
        return np.rad2deg(EA)


def EArate313(t, w, radian_in=False, radian_out=False):
    """
    Calculation of 3-1-3 Euler Angle Rates in a given instant.
    Angular velocities should be inserted in rad/sec.

    radian_in: True if Euler Angles vector are in radian, else false.
    
    radian_out: True if it is desired that the function returns the rates in radians, else false.   
    """
    if not radian_in:
        t = np.deg2rad(t)

    a11 = sin(t[2])
    a12 = cos(t[2])
    a13 = 0
    a21 = cos(t[2])*sin(t[1])
    a22 = -sin(t[2])*sin(t[1])
    a23 = 0
    a31 = -sin(t[2])*cos(t[1])
    a32 = -cos(t[2])*cos(t[1])
    a33 = sin(t[1])

    B = np.array([[a11, a12, a13],
                  [a21, a22, a23],
                  [a31, a32, a33]])
    
    rate = (1/sin(t[1])) * np.matmul(B, w.T)

    if radian_out:
        return rate
    else:
        return np.rad2deg(rate)
