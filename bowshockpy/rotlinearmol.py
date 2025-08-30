"""
This module perform the radiative transfer for a CO rotational transition.
These are the assumptions:

- Negligible population of vibrational excited states
- Negligible centrifugal distortions
- Local Thermodynamical Equilibrium

"""

import numpy as np
import bowshockpy.radtrans as rt
from astropy import constants as const
from astropy import units as u


def gJ(J):
    """
    Degeneracy of the level J at which the measurement was made. For a linear
    molecule, g = 2J + 1

    Parameters
    ----------
    J : int
        Rotational level

    Returns
    -------
    int
        Degeneracy of the level J
    """
    return 2*J + 1


def B0J(J, nu):
    """
    Rigid rotor rotation constant, being nu the frequency for the transition
    J-> J-1. For high J an aditional term is needed

    Parameters
    ----------
    J : int
        Rotational level
    nu : astropy.units.quantity
        Frequency of the transition

    Returns
    -------
    float
        Rigid rotor rotation constant.
    """
    return nu / (2*J)


def EJ(J, B0):
    """
    Energy state of a rigid rotor, neglecting centrifugal distortions

    Parameters
    ----------
    J : int
        Rotational level
    B0 : astropy.units.quantity
        Rotation constant

    Returns
    -------
    astropy.units.quantity
        Energy state of a rotator
    """
    return const.h * B0 * J * (J+1)


def muJ_Jm1(J, mu):
    """
    Computes the dipole moment matrix element squared for rotational transition
    J->J-1

    Parameters
    ----------
    J : int
        Rotational level
    mu : astropy.units.quantity
        Permanent dipole moment of the molecule
    """
    return mu * np.sqrt(J / (2*J + 1))


def tau_linearmol(dNmoldv, J, nu, Tex, mu):
    """
    Computes the opacity as a function of the column density per channel width
    for a rotational transition of a linear molecule

    Parameters
    ----------
    dNmol : astropy.units.quantity
        Column density per velocity bin dv
    J : int
        Rotational level
    nu : astropy.units.quantity
        Frequency
    Tex : astropy.units.quantity
        Excitation temperature
    mu : astropy.units.quantity
        Permanent dipole moment of the molecule

    Returns
    -------
    float
        Opacity
    """
    B0 = B0J(J, nu)
    mu_ul = muJ_Jm1(J, mu)
    tau = rt.tau_func(
        dNmoldv=dNmoldv,
        nu=nu,
        Tex=Tex,
        i=J,
        Ei=EJ,
        gi=gJ,
        mu_ul=mu_ul,
        Ei_args=(B0),
        gi_args=(),
    )
    return tau
