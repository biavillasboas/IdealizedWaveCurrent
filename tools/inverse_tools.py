import numpy as np
from numpy import pi


def az2trig(az):
    """Transforms azimuth to trigonometric angle

    """
    assert ((az<=360) & (az>=0)).all(), "Azimuth out of range"
    theta = np.ma.masked_all(az.shape)
    ind1 = (az>=0.) & (az<=90.)
    ind2 = (az>90.) & (az<=360.)
    theta[ind1] = 90. - az[ind1]
    theta[ind2] =  90. - az[ind2] + 360.
    return theta

def direction_from_to(theta):
    """Transforms direction 'from' to direction 'to'

    """
    direction = np.ma.masked_all(theta.shape)
    ind1 = (theta>=180.)
    ind2 = (theta<180.)
    direction[ind1] = theta[ind1] - 180
    direction[ind2] = theta[ind2] + 180
    return direction

def dir360_dir180(direction):
    """Transforms direction from [0,360] to [-180, 180]

    """
    ind = direction>=180
    direction[ind] = direction[ind]-360
    return direction


def trig2az(theta):
    """Transforms trigonometric angle to azimuth using met convention

    """
    az0 = np.ma.masked_all(theta.shape)
    idx1 = (90>=theta) & (theta>=-90)
    idx2 = (90<theta) & (theta<=180)
    idx3 = (theta<-90) & (theta>=-180)
    az0[idx1] = abs(theta[idx1] - 90)
    az0[idx2] = (90 - theta[idx2]) + 360
    az0[idx3] = abs(theta[idx3]) + 90
    az = az0.copy()
    az[az0<=180] = az0[az0<=180] + 180
    az[az0>180] = az0[az0>180] - 180
    return az

def grad_k(theta, fx, fy):
    """Projects the horizontal gradient to a direction parallel to the wave direction

    """
    return np.cos(theta)*fx + np.sin(theta)*fy

def grad_n(theta, fx, fy):
    """Projects the horizontal gradient to a direction perpendicular to the wave direction

    """
    return - np.sin(theta)*fx + np.cos(theta)*fy
