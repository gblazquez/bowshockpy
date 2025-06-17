import numpy as np

from astropy import units as u
from astropy import constants as const

freq_caract_CO = {'1-0':115.27120180*u.GHz,
                  '2-1':230.53800000*u.GHz,
                  '3-2':345.79598990*u.GHz,
                  '13CO_3-2':330.58796530*u.GHz}

def exp_hnkt(nu, T):
    return np.exp(const.h*nu/(const.k_B*T))

def Bnu_f(nu, T):
    """
    Planck's law. Intensity dimensions (erg s-1 cm-2 Hz-1 sr-1, or Jy sr-1)
    """
    Bnu = 2 * const.h * nu**3 / (const.c**2 * (exp_hnkt(nu,T)-1))
    return Bnu.to(u.Jy) / u.sr

def B0(nu, J):
    """
    Rigid rotor rotation constant, being nu the frequency for the transition
    J-> J-1. For high J an aditional term is needed
    """
    return nu / (2 * J)

def gJ(J):
    """
    Degeneracy of the level J at which the measuremente was made. For a linear
    molecule as CO, g = 2J + 1
    """
    return 2 * J + 1

def Qpart(nu, J, Tex, tol=10**(-15)):
    """
    Calculates the partition function, sum over all states (2J+1)exp(-hBJ(J+1)/kT)
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
    """
    # aj_jm1 = 1.165 * 10**(-11) * mu.to(u.D).value**2 * (J/(2*J+1)) * nu.to(u.GHz).value**3
    # return aj_jm1 * u.s**(-1)
    aj_jm1 = 64 * np.pi**4 * nu.to(u.Hz).value**3 * J * (mu.to(u.D).value*10**(-18))**2 \
             / (3 * const.c.cgs.value**3 * const.h.cgs.value * (2*J+1))
    return aj_jm1 * u.s**(-1)
    #aj_jm1 = 64 * np.pi**4 / (3*const.h*const.c**3) * (mu**2*J/(2*J+1)) * nu**3
    #return aj_jm1

def Ej(nu, J):
    """
    Energy states of a rotator
    """
    return const.h * nu * (J+1) / 2

def coldens(nu, J, mu, Tex, Tbg, Iint):
    """
    Column density for the optically thin case.
    """
    return 8 * np.pi * nu**3 * Qpart(nu,J,Tex) * Iint \
    / (A_j_jm1(nu,J,mu) * gJ(J) * const.c**3 * (exp_hnkt(nu, Tex)-1)
       * np.exp(-Ej(nu,J)/(const.k_B*Tex)) * (Bnu_f(nu,Tex)-Bnu_f(nu,Tbg)))

def Tex(nu21, nu32, I32_21):
    return const.h/(2*const.k_B) * (3*nu21.to(u.Hz)-4*nu32.to(u.Hz)) / np.log(I32_21*(nu21/nu32)**4)

def tau_f(nu, Tex, Tbg, Inu):
    return - np.log(1 - Inu / (Bnu_f(nu,Tex)-Bnu_f(nu,Tbg)))

def Inu_tau(nu, Tex, Tbg, tau):
    return (Bnu_f(nu,Tex)-Bnu_f(nu,Tbg)) * (1 - np.exp(-tau))

def Inu_tau_thin(nu, Tex, Tbg, tau):
    return (Bnu_f(nu,Tex)-Bnu_f(nu,Tbg)) * tau

# def Inu_tau_thin(nu, Tex, Tbg, tau):
#     return (Bnu_f(nu,Tex)-Bnu_f(nu,Tbg)) * tau

def tau_N(nu, J, mu, Tex, dNdv):
    expo = (J+1) * const.h * nu / 2 / const.k_B / Tex
    NJ = (2*J+1) * dNdv / np.exp(expo) / Qpart(nu, J, Tex)
    kk = const.c**3 * A_j_jm1(nu, J, mu) / 8 / np.pi / nu**3
    return kk * NJ * (exp_hnkt(nu, Tex)-1)

def coldens_thick_dv(nu, J, mu, Tex, Tbg, Inu):
    tau = tau_f(nu, Tex, Tbg, Inu)
    return 8 * np.pi * nu**3 * Qpart(nu,J,Tex) * tau \
    / (A_j_jm1(nu,J,mu) * gJ(J) * const.c**3 * (exp_hnkt(nu, Tex)-1)
       * np.exp(-Ej(nu,J)/(const.k_B*Tex)))

def totalmass_thick_dv(nu, J, mu, Tex, Tbg, Inu, distorarea, abund,
                       meanmass):
    if Inu.unit.is_equivalent(u.Jy):
        dist = distorarea
        return dist**2 * coldens_thick_dv(nu, J, mu, Tex, Tbg, Inu)/u.sr \
                * meanmass  / abund
    elif Inu.unit.is_equivalent(u.Jy / u.sr):
        area = distorarea
        return area * coldens_thick_dv(nu, J, mu, Tex, Tbg, Inu) \
                * meanmass  / abund

def coldens_thin_dv(nu, J, mu, Tex, Tbg, Inu):
    """
    Column density for the optically thin case.
    """
    return 8 * np.pi * nu**3 * Qpart(nu,J,Tex) * Inu \
    / (A_j_jm1(nu,J,mu) * gJ(J) * const.c**3 * (exp_hnkt(nu, Tex)-1)
       * np.exp(-Ej(nu,J)/(const.k_B*Tex)) * (Bnu_f(nu,Tex)-Bnu_f(nu,Tbg)))


def totalmass(nu, J, mu, Tex, Tbg, Iint, dist, abund, total_mass=False,
              meanmolweight=2.4):
    """
    From CO line emission, estimates the total mass through the calculation of
    the column density in the optically thin approximation and considering LTE
    conditions.

    Parameters
    ----------
    nu: astropy.units.Unit
        Frequency of the transition
    J: int
        upper rotational state of the transition
    mu: astropy.units.Unit
        dipole moment, should be in Debye units
    Tex: astropy.units.Unit
        Excitation temperature
    Tbg: astropy.units.Unit
        Background temperature
    Iint: astropy.units.Unit
        Integrated flux in u.Jy*u.km/u.s measured over an area, or the
        intensity in u.Jy * u.km / u.s / u.sr. Note that this decision
        determines the quantity of parameter dist. Note also that this quantity
        corresponds to I_nu \delta v, where \delta v is the velocity width of
        the channel or the integration range of a zeroth order moment.
    dist: astropy.units.Unit
        If Iint is u.Jy*u.km/u.s, dist should be the distance to the source; if
        Iint is in u.Jy*u.km/u.s/u.sr, dist should be the physical area of the
        region
    abund: float
        Abundance of CO
    total_mass: bool
        The output is the total mass of H2 (if total_mass is False), or a total
        mass taking into account extra material heavier but less abundant than
        H2 assuming a mean mass of the interestelar medium of meanmolweight
        (default 2.4 times the mass of the hydrogen atom) for molecule (if
        total_mass is True)
    meanmolweight: float
        Mean molecular weight times the mass of a hydrogen atom taking into
        account the abundance of helium and other trace constituents, only
        taken into account if total_mass = True
    """
    h2mass = 2.01588 / (6.023*10**23) * u.g
    meanmass = meanmolweight / (6.023*10**23) * u.g
    if total_mass:
        m = meanmass
    else:
        m = h2mass
    if Iint.unit.is_equivalent(u.Jy * u.km / u.s):
        return dist**2 * coldens(nu, J, mu, Tex, Tbg, Iint)/u.sr * m  / abund
    elif Iint.unit.is_equivalent(u.Jy * u.km / (u.s * u.sr)):
        area = dist
        return area * coldens(nu, J, mu, Tex, Tbg, Iint) * m  / abund

if __name__ == "__main__":
    import argparse

    description = """DESCRIPTION: From CO line emission, estimates the total
                     mass through the calculation of the column density in the
                     optically thin approximation and considering LTE
                     conditions."""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--nu", type=float,
                        help="Frequency of the transition in GHz",
                        default=freq_caract_CO["3-2"].value)
    parser.add_argument("--J", type=int,
                        help="Upper rotational state of the transition",
                        default=3)
    parser.add_argument("--mu", type=float, help="Dipole moment in Debye units",
                        default=0.112)
    parser.add_argument("--Tex", type=float, help="Excitation temperature in K",
                        default=25)
    parser.add_argument("--Tbg", type=float, help="Background temperature in K",
                        default=2.7)
    parser.add_argument("--Iint", type=float,
                        help="Integrated intensity in u.Jy*u.km/u.s",
                        default=1)
    parser.add_argument("--dist", type=float,
                        help="distance to the source in pc",
                        default=300)
    parser.add_argument("--abund", type=float,
                        help="CO abundance",
                        default=10**(-4))
    parser.add_argument("--total_mass", type=str,
                        help="""total mass of H2 (if total mass is false), or a
                        total mass taking into account extra material heavier
                        but less abundant than H2, assuming a mean mass of the
                        interestelar medium equal to the parameter
                        meanmolweight (if total mass is true) """,
                        default="false")
    parser.add_argument("--meanmolweight", type=float,
                        help="""Mean molecular weight times the mass of a
                        hydrogen atom taking into account the abundance of
                        helium and other trace constituents, only taken into
                        account if total_mass = true""",
                        default=2.4)

    args = parser.parse_args()
    total_mass = True if args.total_mass == "true" else False
    totmass = totalmass(nu = args.nu * u.GHz,
                        J = args.J,
                        mu = args.mu * u.D,
                        Tex = args.Tex * u.K,
                        Tbg = args.Tbg * u.K,
                        Iint = args.Iint * u.Jy * u.km / u.s,
                        dist = args.dist * u.pc,
                        abund = args.abund,
                        total_mass = total_mass,
                        meanmolweight = args.meanmolweight).to(u.solMass).value

    print("{:.16f} solar masses".format(totmass))
