"""
This module perform the radiative transfer for a CO rotational transition. This
are the assumptions:

- Negligible population of vibrational excited states
- Negligible centrifugal distortions
- Local Thermodynamical Equilibrium
"""

import numpy as np
from astropy import constants as const
from astropy import units as u

#####
# Linear molecule pure rotational transition, no centrifugal distortion:
#####


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
    nu : astropy.units.quantity
        Frequency of the transition

    Returns
    -------
    float
        Rigid rotor rotation constant.
    """
    return nu / (2 * J)


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
    return const.h * B0 * J * (J + 1)


def muJ_Jm1(J, mu):
    """
    Computes the dipole moment matrix element squared for J->J-1 transition for
    a pure rotational transition.

    Parameters
    ----------
    J : int
        Rotational level
    mu : astropy.units.quantity
        Permanent dipole moment
    """
    return mu * np.sqrt(J / (2 * J + 1))


#####
# General functions
#####


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
        CO abundance relative to molecular hydrogen

    Returns
    -------
    astropy.units.quantity
        CO column density
    """
    return Ntot * abund


def Qpart(Tex, Ei, gi, Ei_args=(), gi_args=(), tol=10 ** (-15)):
    """
    Computes the partition function.

    Parameters
    ----------
    Tex : astropy.units.quntity
        Excitation temperature
    Ei : callable
        Function Ei(i, *Ei_args) to compute the energy at level i
    gi : callable
        Function gi(i, *gi_args) to compute the degeneracy at level i
    Ei_args : tuple, optional
        Extra arguments passed to Ei function
    gi_args : tuple, optional
        Extra arguments passed to gi
    tol : float, optional
        Tolerance at which the summation is stopped, by default 10**(-15)

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
    """
    Computes the column density of a molecule at energy level i.

    Parameters
    ----------
    N : astropy.units.quantity
       Column density of the molecule
    Tex : astropy.units.quantity
        Excitation temperature
    Ei : callable
        Function Ei(i, *Ei_args) to compute the energy at level i
    gi : callable
        Function gi(i, *gi_args) to compute the degeneracy at level i
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
    Nmol_i = g * Nmol / q * np.exp(-e / (const.k_B * Tex))
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
    Acoeff = (
        64
        * np.pi**4
        * nu.cgs.value**3
        / (3 * const.h.cgs.value * const.c.cgs.value**3)
        * (mu_ul.to(u.D).value * 10 ** (-18)) ** 2
    )
    return Acoeff * u.s ** (-1)


def tau_N(
    dNmol,
    dv,
    nu,
    Tex,
    i,
    Ei,
    gi,
    mui,
    Ei_args=(),
    gi_args=(),
    mui_args=(),
):
    """
    Computes the opacity as a function of the column density per channel width

    Parameters
    ----------
    dNmol : astropy.units.quantity
        Column density in a velocity bin dv
    dv : astropy.units.quantity
        Velocity bin width
    nu : astropy.units.quantity
        Frequency
    Tex : astropy.units.quntity
        Excitation temperature
    i : int
        Level
    mui : callable
        Function mui(i, *mui_args) to compute the dipole moment matrix element
        i, i-1
    Ei : callable
        Function Ei(i, *Ei_args) to compute the energy at level i
    gi : callable
        Function gi(i, *gi_args) to compute the degeneracy at level i
    Ei_args : tuple, optional
        Extra arguments passed to Ei function
    gi_args : tuple, optional
        Extra arguments passed to gi
    mui_args : tuple, optional
        Extra arguments passed to mui

    Returns
    -------
    float
        Opacity
    """
    if not isinstance(mui_args, tuple):
        mui_args = (mui_args,)

    dNmol_i = column_density_mol_i(
        Nmol=dNmol, Tex=Tex, i=i, Ei=Ei, gi=gi, Ei_args=Ei_args, gi_args=gi_args
    )

    mm = mui(i, *mui_args)
    A = A_ul(nu, mm)
    tau = const.c**3 * A / (8 * np.pi * nu**3) * (exp_hnkt(nu, Tex) - 1) * dNmol_i / dv
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


def Bnu_f(nu, T):
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


def Inu_tau(nu, Tex, Tbg, tau):
    """
    Computes the intensity through the radiative transfer equation.

    Parameters
    ----------
    nu : astropy.units.quantity
        Frequency of the transition
    Tex : astropy.units.quntity
        Excitation temperature
    Tbg: astropy.units.quantity
        Background temperature
    tau : float
        Opacity

    Returns
    -------
    astropy.units.quantity
        Intensity (energy per unit of area, time, frequency and solid angle)
    """
    return (Bnu_f(nu, Tex) - Bnu_f(nu, Tbg)) * (1 - np.exp(-tau))


def Inu_tau_thin(nu, Tex, Tbg, tau):
    """
    Computes the intensity taking from the radiative transfer equation under the
    optically thin approximation

    Parameters
    ----------
    nu : astropy.units.quantity
        Frequency of the transition
    Tex : astropy.units.quntity
        Excitation temperature
    Tbg: astropy.units.quantity
        Background temperature
    tau : float
        Opacity

    Returns
    -------
    astropy.units.quantity
        Intensity (energy per unit of area, time, frequency and solid angle)
    """
    return (Bnu_f(nu, Tex) - Bnu_f(nu, Tbg)) * tau


# I think the next two function should be deleted.
# Used in next fucntion. Delete also test
def Ntot_opthin_Inudv(nu, J, mu, Tex, Tbg, Inudv):
    """
    Column density for the optically thin case for a given intensity times
    channel velocity width

    Parameters
    ----------
    J : int
        Upper level of the rotational transition
    mu : astropy.units.quantity
        Dipole moment of the molecule
    Tex : astropy.units.quntity
        Excitation temperature
    Tbg: astropy.units.quantity
        Background temperature
    Inudv : astropy.units.quantity
        Intensity (Jy per solid angle units) times channel map width (velocity
        units)

    Returns
    -------
    astropy.units.quantity
        Column density (particles per unit or area)
    """
    Ntot_opthin = (
        8
        * np.pi
        * nu**3
        * Qpart(nu, J, Tex)
        * Inudv
        / (
            A_j_jm1(nu, J, mu)
            * gJ(J)
            * const.c**3
            * (exp_hnkt(nu, Tex) - 1)
            * np.exp(-Ej(nu, J) / (const.k_B * Tex))
            * (Bnu_f(nu, Tex) - Bnu_f(nu, Tbg))
        )
    )
    return Ntot_opthin


# Used in test
def totmass_opthin(nu, J, mu, Tex, Tbg, Inudv, area, meanmolmass, abund):
    """
    Computes the total mass (molecular hydrogen plus heavier components) in the
    assuming optically thin emission.

    Parameters
    ----------
    J : int
        Upper level of the rotational transition
    mu : astropy.units.quantity
        Dipole moment of the molecule
    Tex : astropy.units.quntity
        Excitation temperature
    Tbg : astropy.units.quantity
        Background temperature
    Inudv : astropy.units.quantity
        Intensity (Jy per solid angle units) times channel map width (velocity
        units)
    area : astropy.units.quantity
        Projected area of a pixel
    meanmolmass : float
        Mean molecular mass per hydrogen molecule
    abund : float
        CO abundance relative to molecular hydrogen

    Returns
    -------
    astropy.units.quantity
        Total mass (H2 + heavier elements) in astropy.units.Msun
    """
    nu = nu_J(J)
    Ntot = Ntot_opthin_Inudv(nu, J, mu, Tex, Tbg, Inudv)
    totmass = area * Ntot * meanmolmass * const.m_p / abund
    return totmass.to(u.Msun)
