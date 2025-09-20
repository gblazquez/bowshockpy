"""This module contains the general equations that allows to compute the column
densities and opacities, and to perform the radiative transfer"""

import numpy as np
from astropy import constants as const
from astropy import units as u

import bowshockpy.utils as ut


def column_density_tot(m, area, meanmolmass):
    """
    Computes the total (H2 + heavier components) column density given the mass
    and the projected area

    Parameters
    ----------
    m : float | astropy.units.Quantity
        Mass. If float, units should be Solar Masses.
    area : float | astropy.units.Quantity
        Projected area. If float, the area should be in cm**2.
    meanmolmass : float
        Mean molecular mass per hydrogen molecule

    Returns
    -------
    astropy.units.Quantity
        Total column density (H2 + heavier components)
    """
    m = ut.make_astropy_units(m, u.Msun)
    area = ut.make_astropy_units(area, u.cm**2)
    return m / (meanmolmass * const.m_p * area)


def column_density_mol(Ntot, abund):
    """
    Computes the column density of a molecule given its abundance with respect
    to the H2

    Parameters
    ----------
    Ntot : float | astropy.units.Quantity
        Total column density (H2 + heavier components). If float, the column
        density should be given in particles per cm**2
    abund : float
        Abundance relative to molecular hydrogen

    Returns
    -------
    astropy.units.Quantity
        Column density of the molecule
    """
    Ntot = ut.make_astropy_units(Ntot, u.cm ** (-2))
    return Ntot * abund


def Qpart(Tex, Ei, gi, Ei_args=(), gi_args=(), tol=10 ** (-15)):
    r"""
    Computes the partition function.

    Parameters
    ----------
    Tex : float | astropy.units.quntity
        Excitation temperature. If float, it should be in Kelvin
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
    Tex = ut.make_astropy_units(Tex, u.K)
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
    Nmol : float | astropy.units.Quantity
       Column density of the emitting molecule. If float, should be in cm**(-2)
    Tex : float | astropy.units.Quantity
        Excitation temperature. If float, should be in Kelvin.
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
    astropy.units.Quantity
        Column density of the molecule at energy level i

    """
    Nmol = ut.make_astropy_units(Nmol, u.cm ** (-2))
    Tex = ut.make_astropy_units(Tex, u.K)
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
    nu : float | astropy.units.Quantity
        Frequency of the transition. If float, it should be in GHz.
    mu_ul : float | astropy.units.Quantity
        Dipole moment matrix element u,l. If float, it should be in Debye.
    """

    nu = ut.make_astropy_units(nu, u.GHz)
    mu_ul = ut.make_astropy_units(mu_ul, u.Debye)
    acoeff = (
        64
        * np.pi**4
        * nu.cgs.value**3
        / (3 * const.h.cgs.value * const.c.cgs.value**3)
        * (mu_ul.to(u.Debye).value * 10 ** (-18)) ** 2
    )
    return acoeff * u.s ** (-1)


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
    dNmoldv : float | astropy.units.Quantity
        Column density per velocity bin. If float, it should be given in
        cm**(-2) * s * km**(-1)
    nu : float | astropy.units.Quantity
        Frequency. If float, it should be in GHz.
    Tex : astropy.units.Quantity
        Excitation temperature. If float, it should be in Kelvin.
    i : int
        Level
    mu_ul : float | astropy.units.Quantity
        Dipole moment matrix element i, i-1. If float should be in Debye.
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

    Tex = ut.make_astropy_units(Tex, u.K)
    dNmoldv = ut.make_astropy_units(dNmoldv, u.s / u.km / u.cm**2)
    nu = ut.make_astropy_units(nu, u.GHz)
    mu_ul = ut.make_astropy_units(mu_ul, u.Debye)

    dNmolidv = column_density_mol_i(
        Nmol=dNmoldv,
        Tex=Tex,
        i=i,
        Ei=Ei,
        gi=gi,
        Ei_args=Ei_args,
        gi_args=gi_args,
    )

    A = A_ul(nu, mu_ul)
    tau = (
        const.c**3
        * A
        / (8 * np.pi * nu**3)
        * (exp_hnkt(nu, Tex) - 1)
        * dNmolidv
    )
    return tau


def exp_hnkt(nu, T):
    """
    Computes exp(h nu / k_B/T)

    Parameters
    ----------
    nu : float | astropy.units.Quantity
        Frequency. If float, it should be in GHz

    T : float | astropy.units.Quantity
        Temperature. If float, should be in Kelvin

    Returns
    -------
    float
        exp(h nu / k_B/T)
    """
    T = ut.make_astropy_units(T, u.K)
    nu = ut.make_astropy_units(nu, u.GHz)
    return np.exp(const.h * nu / (const.k_B * T))


def Bnu_func(nu, T):
    """
    Computes the spectral radiance or specific intensity of a Planckian (energy
    per unit of area, time, frequency, and solid angle)

    Parameters
    ----------
    nu : float | astropy.units.Quantity
        Frequency. If float, it should be in GHz.
    T : float | astropy.units.Quantity
        Temperature. If float, it should be in Kelvin.

    Returns
    -------
    Bnu : astropy.units.Quantity
        Spectral radiance in u.Jy / u.sr
    """
    T = ut.make_astropy_units(T, u.K)
    nu = ut.make_astropy_units(nu, u.GHz)
    Bnu = 2 * const.h * nu**3 / (const.c**2 * (exp_hnkt(nu, T) - 1))
    return Bnu.to(u.Jy) / u.sr


def Inu_func(tau, nu, Tex, Tbg):
    """
    Computes the intensity through the radiative transfer equation.

    Parameters
    ----------
    tau : float
        Opacity.
    nu : float | astropy.units.Quantity
        Frequency of the transition. If float, it should be in GHz
    Tex : float | astropy.units.quntity
        Excitation temperature. If float, it should be in Kelvin.
    Tbg: float | astropy.units.Quantity
        Background temperature. If float, it should be in Kelvin

    Returns
    -------
    astropy.units.Quantity
        Intensity (energy per unit of area, time, frequency and solid angle)
    """
    Tex = ut.make_astropy_units(Tex, u.K)
    Tbg = ut.make_astropy_units(Tbg, u.K)
    nu = ut.make_astropy_units(nu, u.GHz)
    return (Bnu_func(nu, Tex) - Bnu_func(nu, Tbg)) * (1 - np.exp(-tau))
