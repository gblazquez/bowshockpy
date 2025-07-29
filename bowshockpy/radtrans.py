import numpy as np

from astropy import units as u
from astropy import constants as const

freq_caract_CO = {
    '1-0': 115.27120180 * u.GHz,
    '2-1': 230.53800000 * u.GHz,
    '3-2': 345.79598990 * u.GHz,
    '13CO_3-2': 330.58796530 * u.GHz
    }

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
    return np.exp(const.h*nu/(const.k_B*T))


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
    Bnu = 2 * const.h * nu**3 / (const.c**2 * (exp_hnkt(nu,T)-1))
    return Bnu.to(u.Jy) / u.sr


def B0(nu, J):
    """
    Rigid rotor rotation constant, being nu the frequency for the transition
    J-> J-1. For high J an aditional term is needed

    Parameters
    ----------
    nu : astropy.units.quantity
        Frequency of the transition
    J : int
        Upper level of the rotational transition

    Returns
    -------
    float
        Rigid rotor rotation constant.
    """
    return nu / (2 * J)


def gJ(J):
    """
    Degeneracy of the level J at which the measurement was made. For a linear molecule as CO, g = 2J + 1

    Parameters
    ----------
    J : int
        Upper level of the rotational transition

    Returns
    -------
    int
        Degeneracy of the level J
    """
    return 2 * J + 1


def Qpart(nu, J, Tex, tol=10**(-15)):
    """
    Computes the partition function, sum over all states (2J+1)exp(-hBJ(J+1)/kT)

    Parameters
    ----------
    nu : astropy.units.quantity
        Frequency of the transition
    J : int
        Upper level of the rotational transition
    Tex : astropy.units.quntity
        Excitation temperature
    tol : float, optional
        Tolerance at which the summation is stopped, by default 10**(-15)

    Returns
    -------
    Qpart : float
        Partition function
    """
    Qs = []
    diff = 1
    j = 0
    while diff > tol:
        q = gJ(j) * np.exp(-const.h*B0(nu,J)*j*(j + 1)/(const.k_B * Tex))
        Qs.append(q)
        diff = np.abs(Qs[-1]-Qs[-2]) if len(Qs)>=2 else 1
        j += 1
    return np.sum(Qs)


def A_j_jm1(nu, J, mu):
    """
    Calculates the spontaneous emission coeffitient for the J -> J-1 transition
    of a molecule with a dipole moment mu.

    Parameters
    ----------
    nu : astropy.units.quantity
        Frequency of the transition
    J : int
        Upper level of the rotational transition
    mu : astropy.units.quantity
        Dipole moment of the molecule

    Returns
    -------
    astropy.units.quantity
        Spontaneous emission coefficient in s**(-1)
    """
    aj_jm1 = 64 * np.pi**4 * nu.to(u.Hz).value**3 * J * (mu.to(u.D).value*10**(-18))**2 \
             / (3 * const.c.cgs.value**3 * const.h.cgs.value * (2*J+1))
    return aj_jm1 * u.s**(-1)


def Ej(nu, J):
    """
    Energy state of a rotator

    Parameters
    ----------
    nu : astropy.units.quantity
        Frequency of the transition
    J : int
        Upper level of the rotational transition
 
    Returns
    -------
    astropy.units.quantity
        Energy state of a rotator
    """
    return const.h * nu * (J+1) / 2


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


def column_density_CO(m, meanmolmass, area, XCO):
    """
    Computes the CO column density given the mass
    and the projected area
    
    Parameters
    ----------
    m : astropy.units.quantity
        Mass
    meanmolmass : float
        Mean molecular mass per hydrogen molecule
    area : astropy.units.quantity
        Projected area
    XCO : float
        CO abundance relative to molecular hydrogen

    Returns
    -------
    astropy.units.quantity
        CO column density
    """
    return column_density_tot(m, meanmolmass, area) * XCO


def tau_N(nu, J, mu, Tex, dNdv):
    """
    Computes the opacity as a function of the column density per channel width
    
    Parameters
    ----------
    nu : astropy.units.quantity
        Frequency of the transition
    J : int
        Upper level of the rotational transition
    Tex : astropy.units.quntity
        Excitation temperature
    mu : astropy.units.quantity
        Dipole moment of the molecule
    dNdv : astropy.units.quantity
        Column density divided by the channel width

    Returns
    -------
    float
        Opacity
    """
    expo = (J+1) * const.h * nu / 2 / const.k_B / Tex
    NJ = (2*J+1) * dNdv / np.exp(expo) / Qpart(nu, J, Tex)
    kk = const.c**3 * A_j_jm1(nu, J, mu) / 8 / np.pi / nu**3
    return kk * NJ * (exp_hnkt(nu, Tex)-1)


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
    return (Bnu_f(nu,Tex)-Bnu_f(nu,Tbg)) * (1 - np.exp(-tau))


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
    return (Bnu_f(nu,Tex)-Bnu_f(nu,Tbg)) * tau


def Ntot_opthin_Inudv(nu, J, mu, Tex, Tbg, Inudv):
    """
    Column density for the optically thin case for a given intensity times
    channel velocity width
    
    Parameters
    ----------
    nu : astropy.units.quantity
        Frequency of the transition
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
    Ntot_opthin =  8 * np.pi * nu**3 * Qpart(nu,J,Tex) * Inudv \
    / (A_j_jm1(nu,J,mu) * gJ(J) * const.c**3 * (exp_hnkt(nu, Tex)-1)
       * np.exp(-Ej(nu,J)/(const.k_B*Tex)) * (Bnu_f(nu,Tex)-Bnu_f(nu,Tbg)))
    return Ntot_opthin

def totmass_opthin(nu, J, mu, Tex, Tbg, Inudv, area, meanmolmass, XCO):
    """
    Computes the total mass (molecular hydrogen plus heavier components) in the
    assuming optically thin emission.
    
    Parameters
    ----------
    nu : astropy.units.quantity
        Frequency of the transition
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
    XCO : float
        CO abundance relative to molecular hydrogen

    Returns
    -------
    astropy.units.quantity
        Total mass (H2 + heavier elements) in astropy.units.Msun
    """
    Ntot = Ntot_opthin_Inudv(nu, J, mu, Tex, Tbg, Inudv)
    totmass = area * Ntot * meanmolmass * const.m_p / XCO 
    return totmass.to(u.Msun)