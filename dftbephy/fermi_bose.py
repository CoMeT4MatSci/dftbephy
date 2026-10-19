import numpy as np


# _FD_XMAX = np.log(np.finfo(float).max)/2
_FD_XMAX = 18.42068 # from BoltzTrap

def fermi(x):
    """
    Fermi–Dirac occupation (unitless)
    -----------------------
        f(x) = 1 / (exp(x) + 1),   x = (eps - mu) / kBT
    """
    x = np.asarray(x, dtype=float)
    f = np.where(x < 0., 1., 0.)
    norm_ind = np.logical_and(x > -_FD_XMAX, x < _FD_XMAX)
    f[norm_ind] = 1/(np.exp(x[norm_ind]) + 1.)
    return f

def dfermi_deps(x):
    x = np.asarray(x, dtype=float)
    f = np.zeros_like(x)
    norm_ind = np.logical_and(x > -_FD_XMAX, x < _FD_XMAX)
    f[norm_ind] = -np.exp(x[norm_ind])/(np.exp(x[norm_ind]) + 1.)**2
    return f

def step(x):
    f = np.where(x < 0., 1., 0.)
    return f