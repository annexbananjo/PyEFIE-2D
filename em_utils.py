import numpy as np
from numpy.lib.scimath import sqrt
from math import pi

""" Define some EM constants
"""
eps0 = 8.8541878176e-12  # F/m
mu0 = 1.25663706144e-06  # H/m
c0 = 299792458  # m/s


def get_lambda(freq, epsr=1.0):
    """ Returns the wavelenght given the frequency (Hz) and relative permittivity
    """
    if freq == 0:
        print("Frequency should not be zero")
        raise ValueError

    lambda0 = c0 / freq
    return lambda0 / sqrt(epsr)


def get_k0(freq):
    """ Returns the wavenumber from frequency k= 2*pi / lambda0
    """
    return (2.0 * pi) / get_lambda(freq)


def get_impedance(epsrc):
    """Given a complex permittivity get the impedance
    """
    return sqrt(mu0/(eps0*epsrc))
