import numpy as np

from scipy.optimize import minimize, minimize_scalar
from scipy.integrate import quad
from scipy.ndimage import rotate, gaussian_filter

import astropy.units as u
from astropy.io import fits

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import cm
from matplotlib import colors

from datetime import datetime

import bowpy.utils as ut
import bowpy.bsutils as bu
import bowpy.comass as comass



class NJ():

    default_kwargs = {
        "rbf_niters": 1000,
    }

    def __init__(self, ps, **kwargs):
        for param in ps:
            setattr(self, param, ps[param])
        for kwarg in self.default_kwargs:
            kwarg_attr = kwargs[kwarg] if kwarg in kwargs else self.default_kwargs[kwarg]
            setattr(self, kwarg, kwarg_attr)

        if self.rbf_obs is None:
            self.rbf = self.rbf_calc()
        else:
            self.rbf = self.rbf_obs
        self.zbf = self.zb_r(self.rbf)

        self.tj = self.zj / self.vj

        self.tj_yr = self.stoyr(self.tj)

        self.rhow = self.rhow_fromintmass_analytical(
            self.rbf,
            self.mass
        )

        self.rhow_gcm3 = self.solMasskm3togcm3(self.rhow)

        self.mp0 = self.mp0_calc(self.rhow)
        self.mp0_solmassyr = self.mp0_calc(self.rhow) / self.stoyr(1)
        self.mpamb_f = self.mpamb_f_calc(self.rhow)
        self.mpamb_f_solmassyr = self.mpamb_f / self.stoyr(1)

    def stoyr(self, value):
        return value * u.s.to(u.yr)

    def solMasskm2togcm2(self, value):
        return value * (u.solMass/u.km**2).to(u.g/u.cm**2)

    def solMasskm3togcm3(self, value):
        return value * (u.solMass/u.km**3).to(u.g/u.cm**3)

    def gamma(self):
        return (self.vj-self.vw) / self.v0

    def rb(self, zb):
        return (self.L0**2*(self.zj-zb))**(1/3)

    def vr(self, zb):
        return self.v0*(1 + 3*self.rb(zb)**2/self.gamma()/self.L0**2)**(-1)

    def vz(self, zb):
        return self.vw + (self.vj-self.vw)*(1+3*self.rb(zb)**2/self.gamma()/self.L0**2)**(-1)

    def vtot(self, zb):
        return np.sqrt(self.vr(zb)**2 + self.vz(zb)**2)

    def alpha(self, zb):
        """
        This is not the alpha of alex apendix!. Alpha from alex apendix is alpha2.
        Note that when vw=0 this angle is constant
        """
        return np.arctan(self.vr(zb) / self.vz(zb))

    def alpha2(self, zb):
        """
        Becareful, in the jupyter notebbok alpha2 depends on rb instead of zb
        """
        return np.arcsin((1+9*self.rb(zb)**4/self.L0**4)**(-0.5))

    def alpha2_rb(self, rb):
        return np.arcsin((1+9*rb**4/self.L0**4)**(-0.5))

    def theta(self, zb):
        return np.arctan(self.rb(zb) / zb)

    def rbf_0(self, rr):
        """
        This is the eq (11) from Tabone et al. (2018), that should be minimize
        to find rbf
        """
        return 1 / self.gamma() * (rr/self.L0)**3 + rr/self.L0 - self.v0*self.zj/self.L0/self.vj

    def rbf_calc(self, ns=None, use_minimize=True):
        if use_minimize:
            bounds = (0, self.rb(0))
            return minimize_scalar(lambda x: np.abs(self.rbf_0(x)),
                                   method="bounded", bounds=bounds).x
        else:
            ns = self.rbf_niters if ns is None else ns
            rrs = np.linspace(0, self.rb(0), ns)
            trials = np.array([np.abs(self.rbf_0(rr)) for rr in rrs])
            return rrs[np.argmin(trials)]

    def zb_r(self, rr):
        return self.zj - rr**3 / self.L0**2

    # def surfdens(self, rr):
    #     return rhow/2/rb * (ps["L0"]**2*gamma(ps)/3 + rb**2)**2 \
    #        / (rb**2*np.cos(alpha2(rb,ps)) + (ps["L0"]**2/3)*np.sin(alpha2(rb,ps)))

    def surfdens(self, zb):
        cosa = np.cos(self.alpha2(zb))
        tana = np.tan(self.alpha2(zb))
        sd = 0.5 * self.rhow * cosa * (self.gamma()*tana+1)**2 * self.rb(zb)
        return sd

    def dr_func(self, zb, dz):
        """
        Differential of r given a differential of z
        """
        return 1/3 * (self.L0 / self.rb(zb))**2 * dz

    def dz_func(self, zb, dr):
        """
        Differential of r given a differential of z
        """
        return 3 * (self.rb(zb) / self.L0)**2 * dr

    def dsurf_func(self, zb, dz, dphi):
        """
        Differential of surface given a differential in z and phi
        """
        # sina = np.sin(self.alpha2(zb))
        sina = (1+9*self.rb(zb)**4/self.L0**4)**(-0.5)
        dr = self.dr_func(zb, dz)
        return self.rb(zb) * dr * dphi / sina

    def dmass_func(self, zb, dz, dphi):
        """
        Differential of mass given a differential in z
        """
        return self.surfdens(zb) * self.dsurf_func(zb, dz, dphi)

    def intmass_analytical(self, rbf):
        uu = rbf / self.L0 * (3/self.gamma())**0.5
        analit_int = uu**5 / 5 + 2*uu**3/3 + uu
        massint = self.rhow * (self.L0/np.sqrt(3))**3 * np.pi * self.gamma()**(5/2) * analit_int
        return massint

    def intmass_numerical(self, r0, rbf, return_residual=False):
            def integrand(rb):
                tana = np.tan(self.alpha2_rb(rb))
                return (self.gamma()*tana+1)**2 / tana * rb**2
            integ = quad(integrand, r0, rbf)
            massint = self.rhow * np.pi * integ[0]
            if return_residual:
                return massint, integ[1]
            else:
                return massint

    def rhow_fromintmass_analytical(self, rb, massint):
        uu = rb / self.L0 * (3/self.gamma())**0.5
        analit_int = uu**5 / 5 + 2*uu**3/3 + uu
        rhow = massint * ((self.L0/np.sqrt(3))**3 * np.pi * self.gamma()**(5/2) * analit_int)**(-1)
        return rhow

    def rhow_fromintmass_sigma_simple(self, R0, Rb, massint, return_residual=False):

        def integrand(rb):
            tana = np.tan(self.alpha2_rb(rb))
            return (self.gamma()*tana+1)**2 / tana * rb**2

        integ = quad(integrand, R0, Rb)
        rhow = massint / np.pi / integ[0]
        if return_residual:
            return rhow, integ[1]
        else:
            return rhow

    def mp0_calc(self, rhow):
        mp0 = rhow * np.pi * self.L0**2 * (self.vj-self.vw)**2 / 3 / self.v0
        return mp0

    def mpamb_f_calc(self, rhow):
        mpamb_f = np.pi * rhow * (self.vj - self.vw) * self.rbf**2
        return mpamb_f



class AJ():

    default_kwargs = {
        "rbf_niters": 1000,
        "rb_niters": 100,
        "zbr_niters": 1000,
        "maxrb": 3,
    }

    def __init__(self, ps, **kwargs):
        for param in ps:
            setattr(self, param.replace("/", "_"), ps[param])
        for kwarg in self.default_kwargs:
            kwarg_attr = kwargs[kwarg] if kwarg in kwargs else self.default_kwargs[kwarg]
            setattr(self, kwarg, kwarg_attr)

        self.zj = (self.xhead - self.xorigin) / np.sin(self.i)
        self.vj = (-self.vhead + self.vsys) / np.cos(self.i)
        self.tj = (self.zj * self.distpc * u.au \
                   / (self.vj * u.km/u.s)).to(u.yr).value
        self.rbf = self.rbf_calc()
        self.zbf = self.zb_r(self.rbf)

    def gamma(self):
        return (self.vj-self.vw) / self.v0

    def rb_0(self, rb, zb):
        return (zb-self.zj) + 1/self.L0**2 * (rb**3 - 3*rb*self.rj**2 + 2*self.rj**3)

    def rb(self, zb, niter=None, return_residual=False, use_minimize=True):
        if use_minimize:
            bounds = (self.rj, self.maxrb)
            return minimize_scalar(lambda x: np.abs(self.rb_0(x,zb)),
                                  method="bounded", bounds=bounds).x
        else:
            niter = niter if niter is not None else self.rb_niters
            rb_trials = np.linspace(self.rj, self.maxrb, niter)
            trials = np.ones(len(rb_trials))
            for trial, rb_trial in enumerate(rb_trials):
                trials[trial] = np.abs(self.rb_0(rb_trial, zb))
            if return_residual:
                return rb_trials[trials.argmin()], trials.min()
            else:
                return rb_trials[trials.argmin()]

    def vr(self, zb):
        return self.v0*(1 + 3*(self.rb(zb)**2-self.rj**2)/self.gamma()/self.L0**2)**(-1)

    def vz(self, zb):
        return self.vw + (self.vj-self.vw)*\
                (1+3*(self.rb(zb)**2-self.rj**2)/self.gamma()/self.L0**2)**(-1)

    def vtot(self, zb):
        return np.sqrt(self.vr(zb)**2 + self.vz(zb)**2)

    def alpha(self, zb):
        """
        This is not the alpha of alex apendix!. Alpha from alex apendix is alpha2.
        Note that when vw=0 this angle is constant
        """
        return np.arctan(self.vr(zb) / self.vz(zb))

    def theta(self, zb):
        return np.arctan(self.rb(zb) / zb)

    def vz_obs(self, zb, phi):
        a = self.alpha(zb)
        return self.vtot(zb) * (np.cos(a)*np.cos(self.i) \
                                - np.sin(a)*np.cos(phi)*np.sin(self.i))

    def x_obs(self, zb, phi):
        return self.rb(zb)*np.cos(phi)*np.cos(self.i) + zb*np.sin(self.i)

    def y_obs(self, zb, phi):
        return self.rb(zb) * np.sin(phi)

    def rbf_0(self, rr):
        """
        This is the eq (A.2) from Tabone et al. (2018), that should be minimize
        to find rbf
        """
        return 1/self.gamma()/self.L0**2 \
            * (rr**3 - self.rj**3 + 3*self.rj**2*(self.rj-rr)) \
            + rr - self.rj - self.zj*self.v0/self.vj

    def rbf_calc(self, ns=None, use_minimize=True):
        if use_minimize:
            bounds = (self.rj, self.rb(0))
            return minimize_scalar(lambda x: np.abs(self.rbf_0(x)), method="bounded",
                                   bounds=bounds).x
        else:
            ns = self.rbf_niters if ns is None else ns
            rrs = np.linspace(self.rj, self.rb(0), ns)
            trials = np.array([np.abs(self.rbf_0(rr)) for rr in rrs])
            return rrs[np.argmin(trials)]

    def zb_0(self, zb, rb):
        return (zb-self.zj) + 1/self.L0**2 * (rb**3 - 3*rb*self.rj**2 + 2*self.rj**3)

    def zb_r(self, rr, niter=None, return_residual=False, use_minimize=True):
        if use_minimize:
            bounds = (0, self.zj)
            return minimize_scalar(lambda x: np.abs(self.zb_0(x,rr)),
                                   method="bounded", bounds=bounds).x
        else:
            niter = niter if niter is not None else self.zbr_niters
            zb_trials = np.linspace(0, self.zj, niter)
            trials = np.ones(len(zb_trials))
            for trial, zb_trial in enumerate(zb_trials):
                trials[trial] = np.abs(self.zb_0(zb_trial, rr))
            if return_residual:
                return zb_trials[trials.argmin()], trials.min()
            else:
                return zb_trials[trials.argmin()]


class ObsModel(NJ):
    def __init__(self, ps, psobs, **kwargs):
        super().__init__(ps, **kwargs)
        for param in psobs:
            setattr(self, param, psobs[param])
        for kwarg in self.default_kwargs:
            kwarg_attr = kwargs[kwarg] if kwarg in kwargs else self.default_kwargs[kwarg]
            setattr(self, kwarg, kwarg_attr)

        self.rj_arcsec = self.km2arcsec(self.rj)
        self.zj_arcsec = self.km2arcsec(self.zj)
        self.L0_arcsec = self.km2arcsec(self.L0)
        self.rbf_arcsec = self.km2arcsec(self.rbf)
        self.zbf_arcsec = self.km2arcsec(self.zbf)

    def km2arcsec(self, value):
        return value * u.km.to(u.au) / self.distpc

    def vzp(self, zb, phi):
        a = self.alpha(zb)
        return self.vtot(zb) * (np.cos(a)*np.cos(self.i) - np.sin(a)*np.cos(phi)*np.sin(self.i))

    def xp(self, zb, phi):
        return self.rb(zb)*np.cos(phi)*np.cos(self.i) + zb*np.sin(self.i)

    def yp(self, zb, phi):
        return self.rb(zb) * np.sin(phi)


class Bowshock2D(ObsModel):
    nzs_init = 200
    nvisos_init = 200
    def __init__(self, ps, psobs, **kwargs):
        super().__init__(ps, psobs, **kwargs)
        self.ps = ps.copy()
        self.psobs = psobs.copy()
        for param in psobs:
            setattr(self, param, psobs[param])
        for kwarg in self.default_kwargs:
            kwarg_attr = kwargs[kwarg] if kwarg in kwargs else self.default_kwargs[kwarg]
            setattr(self, kwarg, kwarg_attr)

        self.nzs = self.nzs_init if "nzs" not in kwargs else kwargs["nzs"]
        self.nvisos = self.nvisos_init if "nvisos" not in kwargs else kwargs["nvisos"]
        self.zsextended = np.array([])
        self.zs = np.array([])
        self.dzs = np.array([])
        self.Rs = np.array([])
        self.vs = np.array([])
        self.vrs = np.array([])
        self.vzs = np.array([])
        self.vs = np.array([])
        self.thetas = np.array([])
        self.xps_phi90 = np.array([])
        self.xps_phi0 = np.array([])
        self.xps_phi180 = np.array([])
        self.vzps_phi0 = np.array([])
        self.vzps_phi90 = np.array([])
        self.vzps_phi180 = np.array([])
        self.maxviso = None
        self.minviso = None
        self.visos = np.array([])
        self.calc_solutions()


    def calc_solutions(self):
        self.zsextended = self.zb_r(
            np.linspace(self.rbf, self.rj, self.nzs)
        )
        self.dzs = (np.diff(self.zsextended[1:]) + np.diff(self.zsextended[:-1])) / 2
        self.zs = self.zsextended[1:-1]

        self.vs = np.array([self.vtot(zb) for zb in self.zs])
        self.Rs = np.array([self.rb(zb) for zb in self.zs])
        self.vrs = np.array([self.vr(zb) for zb in self.zs])
        self.vzs = np.array([self.vz(zb) for zb in self.zs])
        self.vs = np.array([np.sqrt(vrr**2+vzz**2) for vrr, vzz in zip(self.vrs, self.vzs)])
        self.maxvs = np.max(self.vs)
        self.minvs = np.min(self.vs)
        self.thetas = np.array([np.arctan(self.Rs[i] / z)
                       for i, z in enumerate(self.zs)])

        self.xps_phi90 = np.array([self.xp(zb,phi=np.pi/2) for zb in self.zs])
        self.xps_phi0 = np.array([self.xp(zb,phi=0) for zb in self.zs])
        self.xps_phi180 = np.array([self.xp(zb,phi=np.pi) for zb in self.zs])
        self.vzps_phi0 = -np.array([self.vzp(zb,phi=0) for zb in self.zs]) + self.vsys
        self.vzps_phi90 = -np.array([self.vzp(zb,phi=np.pi/2) for zb in self.zs]) + self.vsys
        self.vzps_phi180 = -np.array([self.vzp(zb,phi=np.pi) for zb in self.zs]) + self.vsys

        self.maxviso = np.max(self.vzps_phi180)
        self.minviso = np.min(self.vzps_phi0)
        self.minvzp90 = np.min(self.vzps_phi90)
        self.maxvzp90 = np.max(self.vzps_phi90)
        self.visos = np.linspace(self.maxviso, self.minviso, self.nvisos)

        self.Rs_arcsec = self.km2arcsec(self.Rs)
        self.zs_arcsec = self.km2arcsec(self.zs)
        self.xps_phi0_arcsec = self.km2arcsec(self.xps_phi0)
        self.xps_phi90_arcsec = self.km2arcsec(self.xps_phi90)
        self.xps_phi180_arcsec = self.km2arcsec(self.xps_phi180)

        self.surfdenss = np.array([self.surfdens(zb) for zb in self.zs])
        self.surfdenss_gcm2 = self.solMasskm2togcm2(self.surfdenss)

    def update_params(self, params):
        print("Recalculating model...")
        for param in params:
            if param in self.ps:
                self.ps[param] = params[param]
            elif param in self.psobs:
                self.psobs[param] = params[param]
            setattr(self, param, params[param])
            self.__init__(self.ps, self.psobs)
        print("Updated!\n")



class Bowshock2DPlots(Bowshock2D):
    narrows = 10
    def __init__(self, ps, psobs, **kwargs):
        super().__init__(ps, psobs, **kwargs)
        self.zs_arrows = np.array([])
        self.Rs_arrows = np.array([])
        self.zs_arrows_tip = np.array([])
        self.Rs_arrows_tip = np.array([])
        self.v_arrow_ref = 100
        self.z_arrow_ref_tip = None
        self.R_arrow_ref_tip = None

        self.axs = {}
        self.cbaxs = {}

        self.calc_arrows()
        self.create_axes()
        self.plot()

    def calc_arrows(self):
        idx_arr = int(len(self.zs_arcsec)/self.narrows)
        self.larrow = 1 / np.max([np.max(self.vrs), np.max(self.vzs)])

        self.zs_arrows = self.zs_arcsec[::idx_arr]
        self.Rs_arrows = self.Rs_arcsec[::idx_arr]

        self.zs_arrows_tip = self.zs_arrows + self.larrow * self.vzs[::idx_arr]
        self.Rs_arrows_tip = self.Rs_arrows + self.larrow * self.vrs[::idx_arr]


    def create_axes(self):
        nrow = 2
        ncol = 3
        wspace = 0.4
        hspace = 0.4
        width_ratios = [1] * ncol
        height_ratios = [1] * nrow

        self.fig_model = plt.figure(figsize=(14,8))
        gs = GridSpec(
            nrow, ncol,
            height_ratios=height_ratios,
            width_ratios=width_ratios,
            hspace=hspace, wspace=wspace
        )
        gss = {}
        gss[0] = gs[0, 1].subgridspec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.05,
        )
        gss[1] = gs[0, 2].subgridspec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.05,
        )
        gss[2] = gs[1, 1].subgridspec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.05,
        )
        gss[3] = gs[1, 2].subgridspec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.05,
        )

        self.axs["text"] = plt.subplot(gs[:, 0])
        self.axs[0] = plt.subplot(gss[0][1, 0])
        self.cbaxs[0] = plt.subplot(gss[0][0, 0])
        self.axs[1] = plt.subplot(gss[1][1, 0])
        self.cbaxs[1] = plt.subplot(gss[1][0, 0])
        self.axs[2] = plt.subplot(gss[2][1, 0])
        self.cbaxs[2] = plt.subplot(gss[2][0, 0])
        self.axs[3] = plt.subplot(gss[3][1, 0])
        self.cbaxs[3] = plt.subplot(gss[3][0, 0])
        self.axs["text"].set_axis_off()

    def plot(self):

        showtext = \
            fr"""
            {self.modelname}

            $i = {{{self.i*180/np.pi:.2f}}}^\circ$
            $V_\mathrm{{jet}} = {{{self.vj:.2f}}}$ km/s
            $V_0 = {{{self.v0:.2f}}}$ km/s
            $V_w = {{{self.vw:.2f}}}$ km/s
            $L_0 = {{{self.L0_arcsec:.2f}}}$ arcsec
            $z_\mathrm{{jet}} = {{{self.zj_arcsec:.2f}}}$ arcsec
            $r_\mathrm{{b,f}} = {{{self.rbf_arcsec:.2f}}}$ arcsec
            $t_\mathrm{{jet}} = {{{self.tj_yr:.2f}}}$ yr
            $\rho_w = {{{self.rhow_gcm3*10**20:.2f}}}\times 10^{{-20}}$ g cm$^{{-3}}$
            $\dot{{m}}_0 = {{{self.mp0_solmassyr*10**6:.2f}}}\times10^{{-6}}$ M$_\odot$ yr$^{{-1}}$
            $\dot{{m}}_{{w,f}} = {{{self.mpamb_f_solmassyr*10**6:.2f}}}\times10^{{-6}}$ M$_\odot$ yr$^{{-1}}$
            """

        self.axs["text"].set_axis_off()
        for n, line in enumerate(showtext.split("\n")):
             self.axs["text"].text(0, 0.99-0.07*n, line, fontsize=12,
                              transform=self.axs["text"].transAxes)

        """
        Deprojected shell
        """
        for i, zarcsec in enumerate(self.zs_arcsec):
            c = ut.get_color(
                [self.minvs, self.maxvs],
                self.vs[i],
                "turbo_r",
            )
            self.axs[0].plot(
                zarcsec,
                self.Rs_arcsec[i],
                color=c,
                marker="o",
            )
            self.axs[0].plot(
                zarcsec,
                -self.Rs_arcsec[i],
                color=c,
                marker="o",
            )
        cbar = plt.colorbar(
               cm.ScalarMappable(
                   norm=colors.Normalize(
                       vmax=self.maxvs,
                       vmin=self.minvs),
                   cmap="turbo_r",
               ),
               cax=self.cbaxs[0],
               orientation="horizontal",
        )

        for i in range(len(self.zs_arrows)):
            self.axs[0].annotate(
                "",
                xy=(self.zs_arrows_tip[i], self.Rs_arrows_tip[i]),
                xytext=(self.zs_arrows[i], self.Rs_arrows[i]),
                arrowprops=dict(arrowstyle="->"))
            self.axs[0].annotate(
                "",
                xy=(self.zs_arrows_tip[i], -self.Rs_arrows_tip[i]),
                xytext=(self.zs_arrows[i], -self.Rs_arrows[i]),
                arrowprops=dict(arrowstyle="->"))

        xlims = [np.min(self.zs_arcsec), np.max(self.zs_arcsec)*1.2]
        ylims = [-np.max(self.Rs_arcsec)*1.3, np.max(self.Rs_arcsec)*1.3]
        self.axs[0].set_xlim(xlims)
        self.axs[0].set_ylim(ylims)
        self.z_arrow_ref = xlims[0] + np.diff(xlims)*0.7
        self.R_arrow_ref = ylims[0] + np.diff(ylims)*0.1
        self.z_arrow_ref_tip = self.z_arrow_ref + self.larrow * self.v_arrow_ref
        self.R_arrow_ref_tip = self.R_arrow_ref + self.larrow * 0
        self.axs[0].annotate(
            "",
            xy=(self.z_arrow_ref_tip, self.R_arrow_ref_tip),
            xytext=(self.z_arrow_ref, self.R_arrow_ref),
            arrowprops=dict(arrowstyle="->"))

        self.axs[0].text(
            self.z_arrow_ref+0.075,
            self.R_arrow_ref+0.05,
            f"{self.v_arrow_ref:d} km/s"
        )

        self.axs[0].set_aspect("equal")
        self.axs[0].set_xlabel("$z_b$ [arcsec]")
        self.axs[0].set_ylabel("R [arcsec]")

        self.cbaxs[0].tick_params(
            bottom=False, labelbottom=False,
            top=True, labeltop=True
        )
        self.cbaxs[0].set_xlabel(r"v [km/s]", )
        self.cbaxs[0].xaxis.set_label_position('top')

        """
        Radius vs distance
        """
        for i, xarcsec in enumerate(self.xps_phi90_arcsec):
            c = ut.get_color(
                [self.minvzp90, self.maxvzp90],
                self.vzps_phi90[i],
                "turbo",
            )
            self.axs[1].plot(
                xarcsec,
                self.Rs_arcsec[i],
                color=c,
                marker="o",
            )
            self.axs[1].plot(
                xarcsec,
                -self.Rs_arcsec[i],
                color=c,
                marker="o",
            )
        cbar = plt.colorbar(
               cm.ScalarMappable(
                   norm=colors.Normalize(
                       vmax=self.maxvzp90,
                       vmin=self.minvzp90),
                   cmap="turbo",
               ),
               cax=self.cbaxs[1],
               orientation="horizontal",
        )
        self.axs[1].set_aspect("equal")
        self.axs[1].set_xlabel("$x(\phi=90^\circ)$ [arcsec]")
        self.axs[1].set_ylabel("$y(\phi=90^\circ)$ [arcsec]")

        self.cbaxs[1].tick_params(
            which="both",
            bottom=False, labelbottom=False,
            top=True, labeltop=True
        )
        self.cbaxs[1].set_xlabel(r"$v_{z'}(\phi=90^\circ)$ [km/s]", )
        self.cbaxs[1].xaxis.set_label_position('top')
        self.cbaxs[1].invert_xaxis()

        """
        Surf. Dens.
        """
        minsurfdens, maxsurfdens = np.min(self.surfdenss_gcm2), np.percentile(self.surfdenss_gcm2, 85)
        self.axs[2].plot(
            self.zs_arcsec[self.surfdenss_gcm2<np.percentile(self.surfdenss_gcm2, 85)],
            self.surfdenss_gcm2[self.surfdenss_gcm2<np.percentile(self.surfdenss_gcm2, 85)],
            "b-"
        )
        self.axs[2].set_xlabel(r"$z$ [arcsec]", )
        self.axs[2].set_ylabel(r"$\sigma$ [g cm$^{-2}$]", )

        self.axs[2].set_yscale("log")

        self.cbaxs[2].set_axis_off()


        """
        PV diagram
        """
        for i in range(len(self.zs)):
            c = ut.get_color(
                    [minsurfdens, maxsurfdens],
                    self.surfdenss_gcm2[i],
                    "viridis",
                    norm="log",
            )
            self.axs[3].plot(
                self.xps_phi180_arcsec[i],
                self.vzps_phi180[i],
                marker="o",
                color=c
            )
            self.axs[3].plot(
                self.xps_phi0_arcsec[i],
                self.vzps_phi0[i],
                marker="o",
                color=c
            )
        self.axs[3].invert_yaxis()
        self.axs[3].set_xlabel("$x_p$ [arcsec]")
        self.axs[3].set_ylabel("$v_{z'}$ [km/s]")

        cbar = plt.colorbar(
               cm.ScalarMappable(
                   norm=colors.LogNorm(
                       vmax=maxsurfdens,
                       vmin=minsurfdens),
                   cmap="viridis",
               ),
               cax=self.cbaxs[3],
               orientation="horizontal",
        )

        self.cbaxs[3].tick_params(
            which="both",
            bottom=False, labelbottom=False,
            top=True, labeltop=True
        )
        self.cbaxs[3].set_xlabel(r"$\sigma$ [g cm$^{-2}$]", )
        self.cbaxs[3].xaxis.set_label_position('top')

#
#
        # self.axs[2].set_xlabel("$x_\mathrm{obs}$ (arcsec)")
        # self.axs[2].set_ylabel("$V_\mathrm{LSR}$ (km/s)")
        # self.axs[2].invert_yaxis()
#
        # self.fig_model.canvas.draw_idle()


class BowshockCube(ObsModel):
    def __init__(self, ps, psobs, pscube, **kwargs):
        super().__init__(ps, psobs, **kwargs)
        self.ps = ps.copy()
        self.psobs = psobs.copy()
        for param in pscube:
            setattr(self, param, pscube[param])
        for kwarg in self.default_kwargs:
            kwarg_attr = kwargs[kwarg] if kwarg in kwargs else self.default_kwargs[kwarg]
            setattr(self, kwarg, kwarg_attr)

        self.nrs = None
        self.rs = np.array([])
        self.dr = None
        self.zs = np.array([])
        self.dzs = np.array([])

        self.phis = np.array([])
        self.dphi = None

        self.vs = np.array([])
        self.vchanss = np.array([])

        self.cube = None
        self.cubes = {}


    def cond_populatechan(self, diffv):
        if self.tolfactor_vt is not None:
            return diffv < np.abs(self.vt)*self.tolfactor_vt
        else:
            return True

    def wvzp(self, diffv, dmass):
        normfactor = np.abs(self.chanwidth) / (np.sqrt(np.pi)*np.abs(self.vt))
        em = dmass * np.exp(-(diffv/self.vt)**2) * normfactor
        return em

    def wxpyp(self, chan, vchan, xpix, ypix, dxpix, dypix, vzp, dmass):
        diffv = np.abs(vzp-vchan)
        if self.cond_populatechan(diffv):
            em = self.wvzp(diffv, dmass)
            self.cube[chan, ypix, xpix] += em * (1-dxpix) * (1-dypix)
            self.cube[chan, ypix, xpix+1] += em * dxpix * (1-dypix)
            self.cube[chan, ypix+1, xpix] += em * (1-dxpix) * dypix
            self.cube[chan, ypix+1, xpix+1] += em * dxpix * dypix

    def makecube(self, ):
        if self.verbose:
            ts = []
            print(f"""
 ---------------------
 Starting making cube
 ---------------------
 Channel width: {self.abschanwidth:.3} km/s
 Pixel size: {self.arcsecpix:.4} arcsec/pix\n
 """)

        self.nrs = self.nzs
        self.rs = np.linspace(self.rbf, self.rj, self.nrs)
        self.dr = self.rs[0] - self.rs[1]
        self.zs = self.zb_r(self.rs)
        self.dzs = self.dz_func(self.zb_r(self.rs), self.dr)

        self.phis = np.linspace(0, 2*np.pi, self.nphis+1)[:-1]
        self.dphi = self.phis[1] - self.phis[0]

        self.vs = np.array([self.vtot(zb) for zb in self.zs])
        self.velchans = np.linspace(self.vch0, self.vchf, self.nc)

        self.cube = np.zeros((self.nc, self.nys, self.nxs))
        ci = np.cos(self.i)
        si = np.sin(self.i)


        outsidegrid_warning = True
        for iz, z in enumerate(self.zs):
            if self.verbose:
                t0 = datetime.now()
                print(f"Computing: z = {z:.4f} km, v = {self.vs[iz]:.4} km/s, r = {self.rs[iz]:.4} km")
                print(f"Computing: z = {self.km2arcsec(z):.4f} acsec, v = {self.vs[iz]:.4} km/s, r = {self.km2arcsec(self.rs[iz]):.4} arcsec")

            if iz != len(self.zs)-1:
                dmass = self.dmass_func(z, self.dzs[iz], self.dphi)
            else:
                dmass = self.intmass_analytical(self.dr/2) / self.nphis

            for phi in self.phis:
                xp = self.rs[iz] * np.cos(phi) * ci + z * si
                yp = self.rs[iz] * np.sin(phi)
                vzp = -self.vzp(z, phi)

                xpixcoord = self.km2arcsec(xp) / self.arcsecpix + self.refpix[0]
                ypixcoord = self.km2arcsec(yp) / self.arcsecpix + self.refpix[1]
                xpix = int(xpixcoord)
                ypix = int(ypixcoord)
                if (xpix+1<self.nxs) and (ypix+1<self.nys) and (xpix>0) and (ypix>0):
                    dxpix = xpixcoord - xpix
                    dypix = ypixcoord - ypix

                    for chan, vchan in enumerate(self.velchans):
                        self.wxpyp(chan, vchan, xpix, ypix, dxpix, dypix, vzp, dmass)

                else:
                    if outsidegrid_warning:
                        print("WARNING: Part of the model lie outside the grid!")
                        outsidegrid_warning = False
            if self.verbose:
                tf = datetime.now()
                print(fr"dt = {(tf-t0).total_seconds()*1000:.2f} ms")
                ts.append((tf-t0).total_seconds())
                print(fr"Time = {np.sum(ts):.2f} s")
                print(f"Progress: {iz/len(self.zs)*100:.2f} %")
        if self.verbose:
            print(f"""
------------------------------------
The spectral cube has been generated 
------------------------------------
            """)



class CubeProcessing(BowshockCube):
    default_kwargs = {
    }
    btypes = {
        "m": "mass",
        "I": "Intensity",
        "Ithin": "Intensity",
        "NCO": "CO column density",
        "tau": "Opacity"
    }
    bunits = {
        "m": "SolarMass",
        "I": "Jy/beam",
        "Ithin": "Jy/beam",
        "NCO": "cm-2",
        "tau": "-"
    }
    dos = {
        "s": "add_source",
        "r": "rotate",
        "n": "add_noise",
        "c": "convolve",
    }

    def __init__(self, bscube, mpars, **kwargs):
        self.__dict__ = bscube.__dict__
        for param in mpars:
                   setattr(self, param, mpars[param])
        for kwarg in self.default_kwargs:
            kwarg_attr = kwargs[kwarg] if kwarg in kwargs else self.default_kwargs[kwarg]
            setattr(self, kwarg, kwarg_attr)

        self.cubes = {}
        self.cubes["m"] = self.cube
        self.refpixs = {}
        self.refpixs["m"] = self.refpix
        self.hdrs = {}

        self.areapix_cm = None
        self.beamarea_sr = None
        self.calc_beamarea_sr()
        self.calc_areapix_cm()

    @staticmethod
    def newck(ck, s):
        return f"{ck}_{s}" if "_" not in ck else  ck+s

    @staticmethod
    def q(ck):
        return ck.split("_")[0] if "_" in ck else ck

    def calc_beamarea_sr(self):
        self.beamarea_sr = ut.mb_sa_gaussian_f(
            self.bmaj*u.arcsec,
            self.bmin*u.arcsec
        )

    def calc_areapix_cm(self):
        self.areapix_cm = ((self.arcsecpix * self.distpc * u.au)**2).to(u.cm**2)

    def calc_NCO(self,):
        self.cubes["NCO"] = (
            self.cubes["m"] * u.solMass * self.XCO \
            / self.meanmass / self.areapix_cm
        ).to(u.cm**(-2)).value #* self.NCOfactor
        self.refpixs["NCO"] = self.refpixs["m"]
        if self.verbose:
            print(f"\nCO column densities has been calculated\n")

    def calc_tau(self):
        if "NCO" not in self.cubes:
            self.calc_NCO()
        self.cubes["tau"] = comass.tau_N(
            nu=comass.freq_caract_CO["3-2"],
            J=3, 
            mu=0.112*u.D,
            Tex=self.Tex,
            Tbg=self.Tbg,
            dNdv=self.cubes["NCO"]*u.cm**(-2) / (self.abschanwidth*u.km/u.s),
        ).to("").value
        self.refpixs["tau"] = self.refpixs["m"]
        if self.verbose:
            print(f"\nOpacities has been calculated\n")

    def calc_I(self, opthin=False):
        """
        Intensity in Jy/beam.
        """
        if "tau" not in self.cubes:
            self.calc_tau()
        func_I = comass.Inu_tau_thin if opthin else comass.Inu_tau
        ckI = "Ithin" if opthin else "I"
        self.cubes[ckI] = (func_I(
            nu=comass.freq_caract_CO["3-2"],
            Tex=self.Tex,
            Tbg=self.Tbg,
            tau=self.cubes["tau"],
        )*self.beamarea_sr).to(u.Jy).value
        self.refpixs[ckI] = self.refpixs["m"]
        if self.verbose:
            print(f"\nIntensities has been calculated\n")

    def calc_Ithin(self):
        self.calc_I(self, opthin=True)

    def add_source(self, ck="m", value=None):
        nck = self.newck(ck, "s")
        self.cubes[nck] = np.copy(self.cubes[ck])
        value = value if value is not None else np.max(self.cubes[ck])
        if self.refpixs[ck][1]>=0 and self.refpixs[ck][0]>=0:
            self.cubes[nck][:, self.refpix[1], self.refpix[0]] = value 
        self.refpixs[nck] = self.refpixs[ck]
        if self.verbose:
            print(f"\n{nck} has been added a source in pix [{self.refpixs[nck][0]:.2f}, {self.refpixs[nck][1]:.2f}] pix\n")

    def rotate(self, ck="m", forpv=False):
        nck = self.newck(ck, "r") if ~forpv else self.newck(ck, "R")
        angle = -self.pa-90 if ~forpv else self.pa+90
        self.cubes[nck] = np.zeros_like(self.cubes[ck])
        for chan in range(np.shape(self.cubes[ck])[0]):
            self.cubes[nck][chan] = rotate(
                self.cubes[ck][chan],
                angle=angle,
                reshape=False,
                order=1
            )
        ang = angle * np.pi/180
        centerx = (self.nxs-1)/2
        centery = (self.nys-1)/2
        rp_center_x = self.refpix[0] - centerx
        rp_center_y = self.refpix[1] - centery
        self.refpixs[nck] = [
         +rp_center_x*np.cos(ang) + rp_center_y*np.sin(ang) + centerx,
         -rp_center_x*np.sin(ang) + rp_center_y*np.cos(ang) + centery
        ]
        if self.verbose:
            print(f"\n{nck} has been rotated to a PA = {self.pa} deg\n")
    
    def add_noise(self, ck="m"):
        nck = self.newck(ck, "n")
        self.cubes[nck] = np.zeros_like(self.cubes[ck])
        for chan in range(np.shape(self.cubes[ck])[0]):
            # sigma_noise = self.target_noise * 2 * np.sqrt(np.pi) \
            #          * np.sqrt(self.x_FWHM*self.y_FWHM) / 2.35    
            sigma_noise = np.max(self.cubes[ck]) / self.maxcube2noise
            noise_matrix = np.random.normal(
                0, sigma_noise,
                size=np.shape(self.cubes[ck][chan])
                )
            self.cubes[nck][chan] = self.cubes[ck][chan] + noise_matrix
        self.refpixs[nck] = self.refpixs[ck]

    def convolve(self, ck="m"):
        nck = self.newck(ck, "c")
        self.cubes[nck] = np.zeros_like(self.cubes[ck])
        for chan in range(np.shape(self.cubes[ck])[0]):
            self.cubes[nck][chan] = bu.gaussconvolve(
                self.cubes[ck][chan],
                x_FWHM=self.x_FWHM,
                y_FWHM=self.y_FWHM,
                pa=self.pabeam,
                return_kernel=False,
            )
        self.refpixs[nck] = self.refpixs[ck]
        if self.verbose:
            print(f"\n{nck} has been convolved with a gaussian kernel [{self.x_FWHM:.2f}, {self.y_FWHM:.2f}] pix\n")

    def calc(self, dostrs):
        for ds in dostrs:
            _split = ds.split("_")
            q = _split[0]
            if q not in self.cubes:
                self.__getattribute__(f"calc_{q}")()
            if len(_split) > 1:
                ss = _split[1]
                for i, s in enumerate(ss):
                    ck = q if i==0 else f"{q}_{ss[:i]}"
                    if self.newck(ck, s) not in self.cubes:
                        self.__getattribute__(self.dos[s])(ck=ck)

    def savecube(self, ck):
        self.hdrs[ck] = bu.create_hdr(
            NAXIS1 = np.shape(self.cubes[ck])[0],
            NAXIS2 = np.shape(self.cubes[ck])[1],
            NAXIS3 = np.shape(self.cubes[ck])[2],
            CRVAL1 = self.ra_source_deg,
            CRVAL2 = self.dec_source_deg,
            CRPIX1 = self.refpixs[ck][0],
            CRPIX2 = self.refpixs[ck][1],
            CDELT1 = -self.arcsecpix / 3600,
            CDELT2 = self.arcsecpix  / 3600,
            BTYPE = self.btypes[self.q(ck)],
            BUNIT = self.bunits[self.q(ck)],
            CTYPE3 = "VRAD",
            CRVAL3 = self.velchans[0],
            CDELT3 = self.velchans[1] - self.velchans[0],
            CUNIT3 = "km/s",
            BMAJ = self.bmaj/3600,
            BMIN = self.bmin/3600,
            BPA = self.pabeam
        )
        hdu = fits.PrimaryHDU(self.cubes[ck])
        hdul = fits.HDUList([hdu])
        hdu.header = self.hdrs[ck]
        bu.make_folder(foldername=f'models/{self.modelname}/fits')
        hdul.writeto(f'models/{self.modelname}/fits/{ck}.fits', overwrite=True)                 
        if self.verbose:
            print(f'models/{self.modelname}/fits/{ck}.fits saved')

    def savecubes(self, cks):
        for ck in cks:
            self.savecube(ck)