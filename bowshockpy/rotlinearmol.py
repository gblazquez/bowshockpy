"""This module contains the implementation of the model of a rotational
transition of a linear molecule. The assumptions are: Local Thermodynamical
Equilibrium, negligible population of vibrational excited states, negligible
centrifugal distortions"""

import numpy as np
from astropy import constants as const
from astropy import units as u

import bowshockpy.radtrans as rt
import bowshockpy.utils as ut


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
    return 2 * J + 1


def B0J(J, nu):
    """
    Rigid rotor rotation constant, being nu the frequency for the transition
    J-> J-1. For high J an aditional term is needed

    Parameters
    ----------
    J : int
        Rotational level
    nu : float | astropy.units.Quantity
        Frequency of the transition. If float, it should be in GHz

    Returns
    -------
    float
        Rigid rotor rotation constant.
    """
    nu = ut.make_astropy_units(nu, u.GHz)
    return nu / (2 * J)


def EJ(J, B0):
    """
    Energy state of a rigid rotor, neglecting centrifugal distortions

    Parameters
    ----------
    J : int
        Rotational level
    B0 : float | astropy.units.Quantity
        Rotation constant. If float, it should be in GHz

    Returns
    -------
    astropy.units.Quantity
        Energy state of a rotator
    """
    B0 = ut.make_astropy_units(B0, u.GHz)
    return const.h * B0 * J * (J + 1)


def muJ_Jm1(J, mu):
    """
    Computes the dipole moment matrix element squared for rotational transition
    J->J-1

    Parameters
    ----------
    J : int
        Rotational level
    mu : float | astropy.units.Quantity
        Permanent dipole moment of the molecule. If float, it should be in
        Debye
    """
    mu = ut.make_astropy_units(mu, u.Debye)
    return mu * np.sqrt(J / (2 * J + 1))


def tau_linearmol(dNmoldv, J, nu, Tex, mu):
    """
    Computes the opacity as a function of the column density per channel width
    for a rotational transition of a linear molecule

    Parameters
    ----------
    dNmoldv : float | astropy.units.Quantity
        Column density per velocity bin dv. If float, it should be in s / km /
        cm**2
    J : int
        Rotational level
    nu : float | astropy.units.Quantity
        Frequency. If float, it should be in GHz.
    Tex : float | astropy.units.Quantity
        Excitation temperature. If float, it should be in Kelvin
    mu : float | astropy.units.Quantity
        Permanent dipole moment of the molecule. If float, it should be in
        Debye

    Returns
    -------
    float
        Opacity
    """
    dNmoldv = ut.make_astropy_units(dNmoldv, u.s / u.km / u.cm**2)
    nu = ut.make_astropy_units(nu, u.GHz)
    Tex = ut.make_astropy_units(Tex, u.K)
    mu = ut.make_astropy_units(mu, u.Debye)

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
