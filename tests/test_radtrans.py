import copy

import numpy as np
from astropy import constants as const
from astropy import units as u

from bowshockpy import radtrans as rt

J = 3
nu = 345 * u.GHz
m = 10 ** (-8) * u.Msun
meanmolmass = 2.8
area = 100 * u.au**2
dv = 1 * u.km / u.s
abund = 10 ** (-4)
mu = 0.112 * u.Debye
Tex = 100 * u.K
Nco = rt.column_density_CO(m=m, meanmolmass=meanmolmass, area=area, abund=abund)
dNcodv = (Nco / dv).to(u.cm ** (-2) / u.km * u.s)
Aeinstein_CO10 = 7.45447951542228e-08

tau = rt.tau_N(nu=nu, J=J, mu=mu, Tex=Tex, dNdv=dNcodv).to("")

eq1 = (
    (
        dNcodv
        / tau
        * (mu.to(u.Debye) * 10 ** (-18)) ** 2
        / const.h.cgs
        * J
        / (2 * J + 1)
        * rt.gJ(J=J)
        / rt.Qpart(nu, J, Tex)
        * np.exp(-rt.Ej(nu=nu, J=J) / const.k_B / Tex)
        * (rt.exp_hnkt(nu=nu, T=Tex) - 1)
    )
    .to((u.D) ** 2 / u.erg / u.cm**3)
    .value
)

def test_Aeinstein():
    assert np.isclose(
        rt.A_j_jm1(nu=rt.freq_caract_CO["1-0"], J=1, mu=mu).value, Aeinstein_CO10
    ), "Einstein coefficient A is not well computed"


def test_opacity():
    assert np.isclose(
        eq1, 3 / 8 / np.pi**3
    ), "Opacities are not consistent with the column density"
