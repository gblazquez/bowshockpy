import copy

import numpy as np
from astropy import constants as const
from astropy import units as u

from bowshockpy import radtrans as rt
from bowshockpy import rotlinearmol as rlm

J = 3
nu = 345 * u.GHz
m = 10 ** (-8) * u.Msun
meanmolmass = 2.8
area = 100 * u.au**2
dv = 1 * u.km / u.s
abund = 10 ** (-4)
mu = 0.112 * u.Debye
Tex = 100 * u.K
Ntot = rt.column_density_tot(m=m, area=area, meanmolmass=meanmolmass)
Nco = rt.column_density_mol(Ntot, abund=abund)
dNcodv = (Nco / dv).to(u.cm ** (-2) / u.km * u.s)
Aeinstein_CO10 = 7.45447951542228e-08

freq_caract_CO = {
    "1-0": 115.27120180 * u.GHz,
}

tau = rlm.tau_linearmol(dNmoldv=dNcodv, J=J, nu=nu, Tex=Tex, mu=mu).to("")

eq1 = (
    (
        dNcodv
        / tau
        * (mu.to(u.Debye) * 10 ** (-18)) ** 2
        / const.h.cgs
        * J
        / (2 * J + 1)
        * rlm.gJ(J=J)
        / rt.Qpart(Tex, gi=rlm.gJ, Ei=rlm.EJ, Ei_args=(rlm.B0J(J, nu)))
        * np.exp(-rlm.EJ(J=J, B0=rlm.B0J(J, nu)) / const.k_B / Tex)
        * (rt.exp_hnkt(nu=nu, T=Tex) - 1)
    )
    .to((u.D) ** 2 / u.erg / u.cm**3)
    .value
)


def test_Aeinstein():
    Aeinstein = rt.A_ul(nu=freq_caract_CO["1-0"], mu_ul=rlm.muJ_Jm1(J=1, mu=mu)).value
    assert np.isclose(
        Aeinstein, Aeinstein_CO10
    ), f"Einstein coefficient A is not well computed {Aeinstein}."


def test_opacity():
    assert np.isclose(
        eq1, 3 / 8 / np.pi**3
    ), "Opacities are not consistent with the column density"
