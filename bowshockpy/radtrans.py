"""This module contains the implementation of general equations to compute the
column densities, opacities and to perform the radiative transfer"""

import numpy as np
from astropy import constants as const
from astropy import units as u


def column_density_tot(m, meanmolmass, area):
    """
    Computes the total (H2 + heavier components) column density given the mass
    and the projected area

    Parameters
    ----------
    m : astropy.units.quantity
        Mass
    meanmolmass : astropy.units.quantity
        Mean molecular mass per hydrogen molecule
    area : astropy.units.quantity
        Projected area

    Returns
    -------
    astropy.units.quantity
        Total column density
    """
    return m / (meanmolmass * const.m_p * area)


def column_density_mol(Ntot, abund):
    """
    Computes the column density of a molecule given its abundance with respect
    to the H2

    Parameters
    ----------
    m : astropy.units.quantity
        Mass
    meanmolmass : float
        Mean molecular mass per hydrogen molecule
    area : astropy.units.quantity
        Projected area
    abund : float
        Abundance relative to molecular hydrogen

    Returns
    -------
    astropy.units.quantity
        Column density of the molecule
    """
    return Ntot * abund


def Qpart(Tex, Ei, gi, Ei_args=(), gi_args=(), tol=10 ** (-15)):
    r"""
    Computes the partition function.

    Parameters
    ----------
    Tex : astropy.units.quntity
        Excitation temperature
    Ei : callable
        Function Ei(i, \*Ei_args) to compute the energy at level i
    gi : callable
        Function gi(i, \*gi_args) to compute the degeneracy at level i
    Ei_args : tuple, optional
        Extra arguments passed to Ei function
    gi_args : tuple, optional
        Extra arguments passed to gi
    tol : float, optional
        Tolerance at which the summation is stopped, by default 10\*\*(-15)

    Returns
    -------
    float
        Partition function
    """
    if not isinstance(Ei_args, tuple):
        Ei_args = (Ei_args,)
    if not isinstance(gi_args, tuple):
        gi_args = (gi_args,)

    Qs = []
    diff = tol * 2
    j = 0
    while diff > tol:
        q = gi(j, *gi_args) * np.exp(-Ei(j, *Ei_args) / (const.k_B * Tex))
        Qs.append(q)
        diff = np.abs(Qs[-1] - Qs[-2]) if len(Qs) >= 2 else Qs[-1]
        j += 1
    return np.sum(Qs)


def column_density_mol_i(Nmol, Tex, i, Ei, gi, Ei_args=(), gi_args=()):
    r"""
    Computes the column density of a molecule at energy level i.

    Parameters
    ----------
    N : astropy.units.quantity
       Column density of the molecule
    Tex : astropy.units.quantity
        Excitation temperature
    Ei : callable
        Function Ei(i, \*Ei_args) to compute the energy at level i
    gi : callable
        Function gi(i, \*gi_args) to compute the degeneracy at level i
    Ei_args : tuple, optional
        Extra arguments passed to Ei function
    gi_args : tuple, optional
        Extra arguments passed to gi

    Returns
    -------
    astropy.units.quantity
        Column density of the molecule at energy level i

    """
    if not isinstance(Ei_args, tuple):
        Ei_args = (Ei_args,)
    if not isinstance(gi_args, tuple):
        gi_args = (gi_args,)

    q = Qpart(Tex, Ei=Ei, gi=gi, Ei_args=Ei_args, gi_args=gi_args)
    g = gi(i, *gi_args)
    e = Ei(i, *Ei_args)
    Nmol_i = Nmol * g / q * np.exp(-e / (const.k_B * Tex))
    return Nmol_i


def A_ul(nu, mu_ul):
    """
    Calculates the spontaneous emission coeffitient for the u -> l level
    transition

    Parameters
    ----------
    nu : astropy.units.quantity
        Frequency of the transition
    mu_ul : astropy.units.quantity
        Dipole moment matrix element u,l
    """
    acoeff = (
        64
        * np.pi**4
        * nu.cgs.value**3
        / (3 * const.h.cgs.value * const.c.cgs.value**3)
        * (mu_ul.to(u.D).value * 10 ** (-18)) ** 2
    )
    return acoeff * u.s**(-1)


def tau_func(
    dNmoldv,
    nu,
    Tex,
    i,
    Ei,
    gi,
    mu_ul,
    Ei_args=(),
    gi_args=(),
):
    r"""
    Computes the opacity as a function of the column density per channel width

    Parameters
    ----------
    dNmoldv : astropy.units.quantity
        Column density per velocity bin
    nu : astropy.units.quantity
        Frequency
    Tex : astropy.units.quantity
        Excitation temperature
    i : int
        Level
    mu_ul : astropy.units.quantity
        Dipole moment matrix element i, i-1
    Ei : callable
        Function Ei(i, \*Ei_args) to compute the energy at level i
    gi : callable
        Function gi(i, \*gi_args) to compute the degeneracy at level i
    Ei_args : tuple, optional
        Extra arguments passed to Ei function
    gi_args : tuple, optional
        Extra arguments passed to gi

    Returns
    -------
    float
        Opacity
    """
    dNmolidv = column_density_mol_i(
        Nmol=dNmoldv, Tex=Tex, i=i, Ei=Ei, gi=gi, Ei_args=Ei_args, gi_args=gi_args
    )

    A = A_ul(nu, mu_ul)
    tau = const.c**3 * A / (8 * np.pi * nu**3) * (exp_hnkt(nu, Tex) - 1) * dNmolidv
    return tau


def exp_hnkt(nu, T):
    """
    Computes exp(h nu / k_B/T)

    Parameters
    ----------
    nu : astropy.units.quantity
        Frequency

    T : astropy.units.quantity
        Temperature

    Returns
    -------
    float
        exp(h nu / k_B/T)
    """
    return np.exp(const.h * nu / (const.k_B * T))


def Bnu_func(nu, T):
    """
    Computes the spectral radiance or specific intensity of a Planckian (energy
    per unit of area, time, frequency, and solid angle)

    Parameters
    ----------
    nu : astropy.units.quantity
        Frequency
    T : astropy.units.quantity
        Temperature

    Returns
    -------
    Bnu : astropy.units.quantity
        Spectral radiance in u.Jy / u.sr
    """
    Bnu = 2 * const.h * nu**3 / (const.c**2 * (exp_hnkt(nu, T) - 1))
    return Bnu.to(u.Jy) / u.sr


def Inu_func(tau, nu, Tex, Tbg):
    """
    Computes the intensity through the radiative transfer equation.

    Parameters
    ----------
    tau : float
        Opacity
    nu : astropy.units.quantity
        Frequency of the transition
    Tex : astropy.units.quntity
        Excitation temperature
    Tbg: astropy.units.quantity
        Background temperature

    Returns
    -------
    astropy.units.quantity
        Intensity (energy per unit of area, time, frequency and solid angle)
    """
    return (Bnu_func(nu, Tex) - Bnu_func(nu, Tbg)) * (1 - np.exp(-tau))
