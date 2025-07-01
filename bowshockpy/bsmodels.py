import numpy as np

from scipy.optimize import minimize_scalar
from scipy.integrate import quad
from scipy.ndimage import rotate

import astropy.units as u
from astropy.io import fits

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import cm
from matplotlib import colors

from datetime import datetime

import sys

import bowshockpy.utils as ut
import bowshockpy.comass as comass
import bowshockpy.moments as moments


class NarrowJet():

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
        return (self.vj-self.va) / self.v0

    def rb(self, zb):
        return (self.L0**2*(self.zj-zb))**(1/3)

    def vr(self, zb):
        return self.v0*(1 + 3*self.rb(zb)**2/self.gamma()/self.L0**2)**(-1)

    def vz(self, zb):
        return self.va + (self.vj-self.va)*(1+3*self.rb(zb)**2/self.gamma()/self.L0**2)**(-1)

    def vtot(self, zb):
        return np.sqrt(self.vr(zb)**2 + self.vz(zb)**2)

    def alpha(self, zb):
        """
        This is not the alpha of alex apendix!. Alpha from alex apendix is alpha2.
        Note that when va=0 this angle is constant
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
        mp0 = rhow * np.pi * self.L0**2 * (self.vj-self.va)**2 / 3 / self.v0
        return mp0

    def mpamb_f_calc(self, rhow):
        mpamb_f = np.pi * rhow * (self.vj - self.va) * self.rbf**2
        return mpamb_f


class ObsModel(NarrowJet):
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
        self.vzps_phi0 = -np.array([self.vzp(zb,phi=0) for zb in self.zs])
        self.vzps_phi90 = -np.array([self.vzp(zb,phi=np.pi/2) for zb in self.zs])
        self.vzps_phi180 = -np.array([self.vzp(zb,phi=np.pi) for zb in self.zs])

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
        # nrow = 2
        nrow = 1
        ncol = 3
        wspace = 0.4
        hspace = 0.4
        width_ratios = [1] * ncol
        height_ratios = [1] * nrow

        self.fig_model = plt.figure(figsize=(14,4))
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
        # gss[2] = gs[1, 1].subgridspec(
        #     2, 1,
        #     height_ratios=[0.05, 1],
        #     width_ratios=[1],
        #     hspace=0.05,
        # )
        # gss[3] = gs[1, 2].subgridspec(
        #     2, 1,
        #     height_ratios=[0.05, 1],
        #     width_ratios=[1],
        #     hspace=0.05,
        # )

        self.axs["text"] = plt.subplot(gs[:, 0])
        self.axs[0] = plt.subplot(gss[0][1, 0])
        self.cbaxs[0] = plt.subplot(gss[0][0, 0])
        self.axs[1] = plt.subplot(gss[1][1, 0])
        self.cbaxs[1] = plt.subplot(gss[1][0, 0])
        # self.axs[2] = plt.subplot(gss[2][1, 0])
        # self.cbaxs[2] = plt.subplot(gss[2][0, 0])
        # self.axs[3] = plt.subplot(gss[3][1, 0])
        # self.cbaxs[3] = plt.subplot(gss[3][0, 0])
        self.axs["text"].set_axis_off()

    def plot(self):

        showtext = \
            fr"""
            {self.modelname}

            $i = {{{self.i*180/np.pi:.2f}}}^\circ$
            $v_\mathrm{{vsys}} = {{{self.vsys:.2f}}}$ km/s
            $v_\mathrm{{iws}} = {{{self.vj:.2f}}}$ km/s
            $v_0 = {{{self.v0:.2f}}}$ km/s
            $v_a = {{{self.va:.2f}}}$ km/s
            $L_0 = {{{self.L0_arcsec:.2f}}}$ arcsec
            $z_\mathrm{{iws}} = {{{self.zj_arcsec:.2f}}}$ arcsec
            $r_\mathrm{{b,f}} = {{{self.rbf_arcsec:.2f}}}$ arcsec
            $t_\mathrm{{iws}} = {{{self.tj_yr:.2f}}}$ yr
            $\rho_a = {{{self.rhow_gcm3*10**20:.2f}}}\times 10^{{-20}}$ g cm$^{{-3}}$
            $\dot{{m}}_0 = {{{self.mp0_solmassyr*10**6:.2f}}}\times10^{{-6}}$ M$_\odot$ yr$^{{-1}}$
            $\dot{{m}}_{{a,f}} = {{{self.mpamb_f_solmassyr*10**6:.2f}}}\times10^{{-6}}$ M$_\odot$ yr$^{{-1}}$
            """

        self.axs["text"].set_axis_off()
        for n, line in enumerate(showtext.split("\n")):
             self.axs["text"].text(0, 0.99-0.06*n, line, fontsize=10,
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
        self.axs[0].set_xlabel("$Distance$ [arcsec]")
        self.axs[0].set_ylabel("Radius [arcsec]")

        self.cbaxs[0].tick_params(
            bottom=False, labelbottom=False,
            top=True, labeltop=True
        )
        self.cbaxs[0].set_xlabel(r"velocity [km/s]", )
        self.cbaxs[0].xaxis.set_label_position('top')

        # """
        # Radius vs distance
        # """
        # for i, xarcsec in enumerate(self.xps_phi90_arcsec):
        #     c = ut.get_color(
        #         [self.minvzp90, self.maxvzp90],
        #         self.vzps_phi90[i],
        #         "turbo",
        #     )
        #     self.axs[1].plot(
        #         xarcsec,
        #         self.Rs_arcsec[i],
        #         color=c,
        #         marker="o",
        #     )
        #     self.axs[1].plot(
        #         xarcsec,
        #         -self.Rs_arcsec[i],
        #         color=c,
        #         marker="o",
        #     )
        # cbar = plt.colorbar(
        #        cm.ScalarMappable(
        #            norm=colors.Normalize(
        #                vmax=self.maxvzp90,
        #                vmin=self.minvzp90),
        #            cmap="turbo",
        #        ),
        #        cax=self.cbaxs[1],
        #        orientation="horizontal",
        # )
        # self.axs[1].set_aspect("equal")
        # self.axs[1].set_xlabel("$x(\phi=90^\circ)$ [arcsec]")
        # self.axs[1].set_ylabel("$y(\phi=90^\circ)$ [arcsec]")

        # self.cbaxs[1].tick_params(
        #     which="both",
        #     bottom=False, labelbottom=False,
        #     top=True, labeltop=True
        # )
        # self.cbaxs[1].set_xlabel(r"$v_{z'}(\phi=90^\circ)$ [km/s]", )
        # self.cbaxs[1].xaxis.set_label_position('top')
        # self.cbaxs[1].invert_xaxis()

        # """
        # Surf. Dens.
        # """
        # minsurfdens, maxsurfdens = np.min(self.surfdenss_gcm2), np.percentile(self.surfdenss_gcm2, 85)
        # self.axs[2].plot(
        #     self.zs_arcsec[self.surfdenss_gcm2<np.percentile(self.surfdenss_gcm2, 85)],
        #     self.surfdenss_gcm2[self.surfdenss_gcm2<np.percentile(self.surfdenss_gcm2, 85)],
        #     "b-"
        # )
        # self.axs[2].set_xlabel(r"$z$ [arcsec]", )
        # self.axs[2].set_ylabel(r"$\sigma$ [g cm$^{-2}$]", )

        # self.axs[2].set_yscale("log")

        # self.cbaxs[2].set_axis_off()


        """
        PV diagram
        """

        minsurfdens = np.min(self.surfdenss_gcm2)
        maxsurfdens = np.percentile(self.surfdenss_gcm2, 85)
        for i in range(len(self.zs)):
            c = ut.get_color(
                    [minsurfdens, maxsurfdens],
                    self.surfdenss_gcm2[i],
                    "viridis",
                    norm="log",
            )
            self.axs[1].plot(
                self.xps_phi180_arcsec[i],
                self.vzps_phi180[i] + self.vsys,
                marker="o",
                color=c
            )
            self.axs[1].plot(
                self.xps_phi0_arcsec[i],
                self.vzps_phi0[i] + self.vsys,
                marker="o",
                color=c
            )
        allvelsarray = np.array([self.vzps_phi0, self.vzps_phi180]).ravel()
        argmaxvelpv = np.argmax(np.abs(allvelsarray))
        if allvelsarray[argmaxvelpv]<0:
            self.axs[1].invert_yaxis()
        else:
            pass
        self.axs[1].set_xlabel("Projected distance [arcsec]")
        self.axs[1].set_ylabel("LOS velocity [km/s]")

        cbar = plt.colorbar(
               cm.ScalarMappable(
                   norm=colors.LogNorm(
                       vmax=maxsurfdens,
                       vmin=minsurfdens),
                   cmap="viridis",
               ),
               cax=self.cbaxs[1],
               orientation="horizontal",
        )

        self.cbaxs[1].tick_params(
            which="both",
            bottom=False, labelbottom=False,
            top=True, labeltop=True
        )
        self.cbaxs[1].set_xlabel(r"Surface density [g cm$^{-2}$]", )
        self.cbaxs[1].xaxis.set_label_position('top')

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
        self.velchans = np.array([])

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

    def makecube(self, fromcube=None):
        if self.verbose:
            ts = []
            print("\nComputing masses...")

        self.nrs = self.nzs
        self.rs = np.linspace(self.rbf, self.rj, self.nrs)
        self.dr = self.rs[0] - self.rs[1]
        self.zs = self.zb_r(self.rs)
        self.dzs = self.dz_func(self.zb_r(self.rs), self.dr)

        self.phis = np.linspace(0, 2*np.pi, self.nphis+1)[:-1]
        self.dphi = self.phis[1] - self.phis[0]

        self.vs = np.array([self.vtot(zb) for zb in self.zs])
        self.velchans = np.linspace(self.vch0, self.vchf, self.nc)

        if fromcube is None:
            self.cube = np.zeros((self.nc, self.nys, self.nxs))
        elif (fromcube is not None) and np.shape(fromcube)==((self.nc, self.nys, self.nxs)):
            self.cube = fromcube
        else:
            sys.exit(f"""
The provided cube into which the model is to be build has dimensions
{np.shape(fromcube)} but the dimensions of the desired model cube is {(self.nc,
self.nys, self.nxs)}. Please, provide a cube with the right dimensions or do not
provide any cube.
""")

        ci = np.cos(self.i)
        si = np.sin(self.i)
        cpa = np.cos(self.pa)
        spa = np.sin(self.pa)

        outsidegrid_warning = True
        ut.progressbar_bowshock(0, self.nzs, length=50, timelapsed=0, intervaltime=0)
        for iz, z in enumerate(self.zs):
            if self.verbose:
                t0 = datetime.now()
#                print(f"Computing: z = {z:.4f} km, v = {self.vs[iz]:.4} km/s, r = {self.rs[iz]:.4} km")
#                print(f"Computing: z = {self.km2arcsec(z):.4f} acsec, v = {self.vs[iz]:.4} km/s, r = {self.km2arcsec(self.rs[iz]):.4} arcsec")

            if iz != len(self.zs)-1:
                dmass = self.dmass_func(z, self.dzs[iz], self.dphi)
            else:
                dmass = self.intmass_analytical(self.dr/2) / self.nphis

            for phi in self.phis:
                _xp = self.rs[iz] * np.sin(phi)
                _yp = self.rs[iz] * np.cos(phi) * ci + z * si
                xp = _xp * cpa - _yp * spa
                yp = _xp * spa + _yp * cpa
                vzp = -self.vzp(z, phi)
                vlsr = vzp + self.vsys

                xpixcoord = self.km2arcsec(xp) / self.arcsecpix + self.refpix[0]
                ypixcoord = self.km2arcsec(yp) / self.arcsecpix + self.refpix[1]
                xpix = int(xpixcoord)
                ypix = int(ypixcoord)
                if (xpix+1<self.nxs) and (ypix+1<self.nys) and (xpix>0) and (ypix>0):
                    dxpix = xpixcoord - xpix
                    dypix = ypixcoord - ypix

                    for chan, vchan in enumerate(self.velchans):
                        self.wxpyp(chan, vchan, xpix, ypix, dxpix, dypix, vlsr, dmass)

                else:
                    if outsidegrid_warning:
                        print("WARNING: Part of the model lie outside the grid! Try to make the image bigger by increasing the number of pixels (parameters nxs and nys), or try to change the reference pixel where the physical center (the source) is found (parameter refpix)\n")
                        outsidegrid_warning = False
            if self.verbose:
                tf = datetime.now()
                intervaltime = (tf-t0).total_seconds()
                ts.append(intervaltime)
                ut.progressbar_bowshock(iz+1, self.nzs, np.sum(ts), intervaltime, length=50)
                # print(f"\n dt = {(tf-t0).total_seconds()*1000:.2f} ms")
                # ts.append((tf-t0).total_seconds())
                # print(f"Time = {np.sum(ts):.2f} s")
                # print(f"Progress: {iz/len(self.zs)*100:.2f} %")



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
    momtol_clipping = 10**(-4)

    def __init__(self, bscube, mpars, **kwargs):
        self.__dict__ = bscube.__dict__
        for param in mpars:
                   setattr(self, param, mpars[param])
        for kwarg in self.default_kwargs:
            kwarg_attr = kwargs[kwarg] if kwarg in kwargs else self.default_kwargs[kwarg]
            setattr(self, kwarg, kwarg_attr)

        self.cubes = {}
        self.cubes["m"] = self.cube
        self.sigma_noises = {}
        self.sigma_noises["m"] = 0
        self.noisychans = {}
        self.noisychans["m"] = np.zeros_like(self.cube[0])
        self.refpixs = {}
        self.refpixs["m"] = self.refpix
        self.hdrs = {}
        self.listmompvs = []

        self.areapix_cm = None
        self.beamarea_sr = None
        self.calc_beamarea_sr()
        self.calc_areapix_cm()

    @staticmethod
    def newck(ck, s):
        return f"{ck}_{s}" if "_" not in ck else ck+s

    @staticmethod
    def q(ck):
        return ck.split("_")[0] if "_" in ck else ck

    def getunitlabel(self, ck):
        return f"{self.btypes[self.q(ck)]} [{self.bunits[self.q(ck)]}]"

    def calc_beamarea_sr(self):
        self.beamarea_sr = ut.mb_sa_gaussian_f(
            self.bmaj*u.arcsec,
            self.bmin*u.arcsec
        )

    def calc_areapix_cm(self):
        self.areapix_cm = ((self.arcsecpix * self.distpc * u.au)**2).to(u.cm**2)

    def calc_NCO(self,):
        if self.verbose:
            print(f"\nComputing column densities...")
        self.cubes["NCO"] = (
            self.cubes["m"] * u.solMass * self.XCO \
            / self.meanmass / self.areapix_cm
        ).to(u.cm**(-2)).value #* self.NCOfactor
        self.refpixs["NCO"] = self.refpixs["m"]
        self.noisychans["NCO"] = self.noisychans["m"]
        self.sigma_noises["NCO"] = self.sigma_noises["m"]
        if self.verbose:
            print(f"CO column densities has been calculated\n")

    def calc_tau(self):
        if "NCO" not in self.cubes:
            self.calc_NCO()
        if self.verbose:
            print(f"\nComputing opacities...")
        self.cubes["tau"] = comass.tau_N(
            nu=comass.freq_caract_CO[self.J],
            J=float(self.J[0]),
            mu=0.112*u.D,
            Tex=self.Tex,
            dNdv=self.cubes["NCO"]*u.cm**(-2) / (self.abschanwidth*u.km/u.s),
        ).to("").value
        self.refpixs["tau"] = self.refpixs["m"]
        self.noisychans["tau"] = self.noisychans["m"]
        self.sigma_noises["tau"] = self.sigma_noises["m"]
        if self.verbose:
            print(f"Opacities has been calculated\n")

    def calc_I(self, opthin=False):
        """
        Intensity in Jy/beam.
        """
        if "tau" not in self.cubes:
            self.calc_tau()
        if self.verbose:
            print(f"\nComputing intensities...")
        func_I = comass.Inu_tau_thin if opthin else comass.Inu_tau
        ckI = "Ithin" if opthin else "I"
        self.cubes[ckI] = (func_I(
            nu=comass.freq_caract_CO[self.J],
            Tex=self.Tex,
            Tbg=self.Tbg,
            tau=self.cubes["tau"],
        )*self.beamarea_sr).to(u.Jy).value
        self.refpixs[ckI] = self.refpixs["m"]
        self.noisychans[ckI] = self.noisychans["m"]
        self.sigma_noises[ckI] = self.sigma_noises["m"]
        if self.verbose:
            print(f"Intensities has been calculated\n")

    def calc_Ithin(self):
        self.calc_I(opthin=True)

    def add_source(self, ck="m", value=None):
        nck = self.newck(ck, "s")
        if self.verbose:
            print(f"\nAdding source to {nck}...")
        self.cubes[nck] = np.copy(self.cubes[ck])
        value = value if value is not None else np.max(self.cubes[ck])
        if self.refpixs[ck][1]>=0 and self.refpixs[ck][0]>=0:
            self.cubes[nck][:, self.refpix[1], self.refpix[0]] = value
        self.refpixs[nck] = self.refpixs[ck]
        self.noisychans[nck] = self.noisychans[ck]
        self.sigma_noises[nck] = self.sigma_noises[ck]
        if self.verbose:
            print(f"{nck} has been added a source in pix [{self.refpixs[nck][0]:.2f}, {self.refpixs[nck][1]:.2f}] pix\n")

    def rotate(self, ck="m", forpv=False):
        nck = self.newck(ck, "r") if not forpv else self.newck(ck, "R")
        if self.verbose:
            if forpv:
                print(f"\nRotatng {nck} in order to compute the PV diagram...")
            else:
                print(f"\nRotatng {nck}...")
        # before allowing rotation of the model and not the cube
        # angle = -self.pa-90 if not forpv else self.pa+90
        # after allowing the model to be rotated
        angle = -self.parot if not forpv else self.papv + 90
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
        rp_center_x = self.refpixs[ck][0] - centerx
        rp_center_y = self.refpixs[ck][1] - centery
        self.refpixs[nck] = [
         +rp_center_x*np.cos(ang) + rp_center_y*np.sin(ang) + centerx,
         -rp_center_x*np.sin(ang) + rp_center_y*np.cos(ang) + centery
        ]
        self.noisychans[nck] = rotate(
            self.noisychans[ck],
            angle=angle,
            reshape=False,
            order=1,
        )
        self.sigma_noises[nck] = self.sigma_noises[ck]
        if self.verbose:
            print(f"{nck} has been rotated {angle} deg to compute the PV diagram\n")

    def add_noise(self, ck="m"):
        nck = self.newck(ck, "n")
        if self.verbose:
            print(f"\nAdding noise to {nck}...")
        self.cubes[nck] = np.zeros_like(self.cubes[ck])
        for chan in range(np.shape(self.cubes[ck])[0]):
            # sigma_noise = self.target_noise * 2 * np.sqrt(np.pi) \
            #          * np.sqrt(self.x_FWHM*self.y_FWHM) / 2.35
            sigma_noise = self.sigma_beforeconv if self.sigma_beforeconv is not None\
                else np.max(self.cubes[ck]) * self.maxcube2noise
            noise_matrix = np.random.normal(
                0, sigma_noise,
                size=np.shape(self.cubes[ck][chan])
                )
            self.cubes[nck][chan] = self.cubes[ck][chan] + noise_matrix
        self.refpixs[nck] = self.refpixs[ck]
        self.noisychans[nck] = noise_matrix
        self.sigma_noises[nck] = sigma_noise
        if self.verbose:
            print(f"Noise added to {nck}\n")

    def convolve(self, ck="m"):
        nck = self.newck(ck, "c")
        if self.verbose:
            print(f"\nConvolving {nck}... ")
        self.cubes[nck] = np.zeros_like(self.cubes[ck])
        for chan in range(np.shape(self.cubes[ck])[0]):
            self.cubes[nck][chan] = ut.gaussconvolve(
                self.cubes[ck][chan],
                x_FWHM=self.x_FWHM,
                y_FWHM=self.y_FWHM,
                pa=self.pabeam,
                return_kernel=False,
            )
        self.refpixs[nck] = self.refpixs[ck]
        self.noisychans[nck] = ut.gaussconvolve(
            self.noisychans[ck],
            x_FWHM=self.x_FWHM,
            y_FWHM=self.y_FWHM,
            pa=self.pabeam,
            return_kernel=False,
        )
        self.sigma_noises[nck] = np.std(self.noisychans[nck])
        if self.verbose:
            print(f"{nck} has been convolved with a gaussian kernel [{self.x_FWHM:.2f}, {self.y_FWHM:.2f}] pix and a PA of {self.pabeam:.2f}deg\n")
            if "n" in nck:
                print(f"The rms of the convolved image is {self.sigma_noises[nck]:.5} {self.bunits[self.q(nck)]}\n")

    def _useroutputcube2dostr(self, userdic):
        dictrad = {
            "mass": "m",
            "intensity": "I",
            "intensity_opthin": "Ithin",
            "CO_column_density": "NCO",
            "opacity": "tau",
            "add_source": "s",
            "rotate": "r",
            "add_noise": "n",
            "convolve": "c",
        }
        dostrs = []
        for userkey in userdic:
            q = dictrad[userkey]
            ops = userdic[userkey]
            calcmompv = "moments_and_pv" in ops
            if calcmompv:
                if len(ops)>1:
                    ss = "".join([dictrad[s_user] for s_user in
                                  userdic[userkey] if s_user!="moments_and_pv"])
                    dostr = [f"{q}_{ss}"]
                if len(ops)==1:
                    dostr = [f"{q}"]
                dostrs += dostr
                self.listmompvs += dostr
            else:
                if len(ops)!=0:
                    ss = "".join([dictrad[s_user] for s_user in userdic[userkey]])
                    dostrs += [f"{q}_{ss}"]
                else:
                    dostrs += [f"{q}"]
        return dostrs

    def calc(self, userdic):
        dostrs = self._useroutputcube2dostr(userdic)
        print(dostrs)
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
        self.hdrs[ck] = ut.create_hdr(
            NAXIS1 = np.shape(self.cubes[ck])[0],
            NAXIS2 = np.shape(self.cubes[ck])[1],
            NAXIS3 = np.shape(self.cubes[ck])[2],
            CRVAL1 = self.ra_source_deg,
            CRVAL2 = self.dec_source_deg,
            CRPIX1 = self.refpixs[ck][0] + 1,
            CRPIX2 = self.refpixs[ck][1] + 1,
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
        ut.make_folder(foldername=f'models/{self.modelname}/fits')
        hdul.writeto(f'models/{self.modelname}/fits/{ck}.fits', overwrite=True)
        if self.verbose:
            print(f'models/{self.modelname}/fits/{ck}.fits saved')

    def savecubes(self, userdic):
        cks = self._useroutputcube2dostr(userdic)
        for ck in cks:
            self.savecube(ck)

    def pvalongz(self, ck, halfwidth=0, save=False):
        pvimage = moments.pv(
            self.cubes[ck],
            int(self.refpixs[ck][1]),
            halfwidth=halfwidth,
            axis=1
        )
        if save:
            hdrpv = ut.create_hdr(
                NAXIS1 = np.shape(self.cubes[ck])[1],
                NAXIS2 = np.shape(self.cubes[ck])[0],
                BMAJ = self.bmaj/3600,
                BMIN = self.bmin/3600,
                BPA = self.pabeam,
                BTYPE = self.btypes[self.q(ck)],
                BUNIT = self.bunits[self.q(ck)],
                CTYPE1 = "OFFSET",
                CRVAL1 = 0,
                CDELT1 = self.arcsecpix,
                CRPIX1 = self.refpixs[ck][0] + 1,
                CUNIT1 = "arcsec",
                CTYPE2 = "VELOCITY",
                CRVAL2 = self.velchans[0],
                CDELT2 = self.chanwidth,
                CRPIX2 = 1,
                CUNIT2 = "km/s"
            )
            hdu = fits.PrimaryHDU(pvimage)
            hdul = fits.HDUList([hdu])
            hdu.header = hdrpv
            ut.make_folder(foldername=f'models/{self.modelname}/fits')
            hdul.writeto(f'models/{self.modelname}/fits/{ck}_pv.fits', overwrite=True)
            if self.verbose:
                print(f'models/{self.modelname}/fits/{ck}_pv.fits saved')
        return pvimage

    def sumint(self, ck, chan_range=None, save=False):
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        sumintimage = moments.sumint(
            self.cubes[ck],
            chan_range=chan_range
        )
        if save:
            hdr = ut.create_hdr(
                NAXIS1 = np.shape(self.cubes[ck])[1],
                NAXIS2 = np.shape(self.cubes[ck])[0],
                BMAJ = self.bmaj/3600,
                BMIN = self.bmin/3600,
                BPA = self.pabeam,
                BTYPE = self.btypes[self.q(ck)],
                BUNIT = self.bunits[self.q(ck)],
                CTYPE1 = "OFFSET",
                CRVAL1 = 0,
                CDELT1 = self.arcsecpix,
                CRPIX1 = self.refpixs[ck][0] + 1,
                CUNIT1 = "arcsec",
                CTYPE2 = "OFFSET",
                CRVAL2 = 0,
                CDELT2 = self.arcsecpix,
                CRPIX2 = self.refpixs[ck][1] + 1,
                CUNIT2 = "arcsec"
            )
            hdu = fits.PrimaryHDU(sumintimage)
            hdul = fits.HDUList([hdu])
            hdu.header = hdr
            ut.make_folder(foldername=f'models/{self.modelname}/fits')
            hdul.writeto(f'models/{self.modelname}/fits/{ck}_sumint.fits', overwrite=True)
            if self.verbose:
                print(f'models/{self.modelname}/fits/{ck}_sumint.fits saved')
        return sumintimage

    def mom0(self, ck, chan_range=None, save=False):
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        chan_vels = self.velchans[chan_range[0]:chan_range[-1]]
        mom0 = moments.mom0(
            self.cubes[ck],
            chan_vels=chan_vels,
            chan_range=chan_range,
            )
        if save:
            hdr = ut.create_hdr(
                NAXIS1 = np.shape(self.cubes[ck])[1],
                NAXIS2 = np.shape(self.cubes[ck])[0],
                BMAJ = self.bmaj/3600,
                BMIN = self.bmin/3600,
                BPA = self.pabeam,
                BTYPE = self.btypes[self.q(ck)],
                BUNIT = self.bunits[self.q(ck)],
                CTYPE1 = "OFFSET",
                CRVAL1 = 0,
                CDELT1 = self.arcsecpix,
                CRPIX1 = self.refpixs[ck][0] + 1,
                CUNIT1 = "arcsec",
                CTYPE2 = "OFFSET",
                CRVAL2 = 0,
                CDELT2 = self.arcsecpix,
                CRPIX2 = self.refpixs[ck][1] + 1,
                CUNIT2 = "arcsec"
            )
            hdu = fits.PrimaryHDU(mom0)
            hdul = fits.HDUList([hdu])
            hdu.header = hdr
            ut.make_folder(foldername=f'models/{self.modelname}/fits')
            hdul.writeto(f'models/{self.modelname}/fits/{ck}_mom0.fits', overwrite=True)
            if self.verbose:
                print(f'models/{self.modelname}/fits/{ck}_mom0.fits saved')
        return mom0

    def mom1(self, ck, chan_range=None, clipping=0, save=False):
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        chan_vels = self.velchans[chan_range[0]:chan_range[-1]]
        cube_clipped = np.copy(self.cubes[ck])
        clipping = clipping if clipping != 0 \
            else self.momtol_clipping * np.max(self.cubes[ck])
        cube_clipped[cube_clipped<clipping] = 0
        mom1 = np.nan_to_num(
                moments.mom1(
                    cube_clipped,
                    chan_vels=chan_vels,
                    chan_range=chan_range
                )
            )
        if save:
            hdr = ut.create_hdr(
                NAXIS1 = np.shape(self.cubes[ck])[1],
                NAXIS2 = np.shape(self.cubes[ck])[0],
                BMAJ = self.bmaj/3600,
                BMIN = self.bmin/3600,
                BPA = self.pabeam,
                BTYPE = self.btypes[self.q(ck)],
                BUNIT = self.bunits[self.q(ck)],
                CTYPE1 = "OFFSET",
                CRVAL1 = 0,
                CDELT1 = self.arcsecpix,
                CRPIX1 = self.refpixs[ck][0] + 1,
                CUNIT1 = "arcsec",
                CTYPE2 = "OFFSET",
                CRVAL2 = 0,
                CDELT2 = self.arcsecpix,
                CRPIX2 = self.refpixs[ck][1] + 1,
                CUNIT2 = "arcsec"
            )
            hdu = fits.PrimaryHDU(mom1)
            hdul = fits.HDUList([hdu])
            hdu.header = hdr
            ut.make_folder(foldername=f'models/{self.modelname}/fits')
            hdul.writeto(f'models/{self.modelname}/fits/{ck}_mom1.fits', overwrite=True)
            if self.verbose:
                print(f'models/{self.modelname}/fits/{ck}_mom1.fits saved')
        return mom1

    def mom2(self, ck, chan_range=None, clipping=0, save=False):
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        chan_vels = self.velchans[chan_range[0]:chan_range[-1]]
        cube_clipped = np.copy(self.cubes[ck])
        clipping = clipping if clipping != 0 \
            else self.momtol_clipping * np.max(self.cubes[ck])
        cube_clipped[cube_clipped<clipping] = 0
        mom2 = np.nan_to_num(
                moments.mom2(
                    cube_clipped,
                    chan_vels=chan_vels,
                    chan_range=chan_range
                )
            )
        if save:
            hdr = ut.create_hdr(
                NAXIS1 = np.shape(self.cubes[ck])[1],
                NAXIS2 = np.shape(self.cubes[ck])[0],
                BMAJ = self.bmaj/3600,
                BMIN = self.bmin/3600,
                BPA = self.pabeam,
                BTYPE = self.btypes[self.q(ck)],
                BUNIT = self.bunits[self.q(ck)],
                CTYPE1 = "OFFSET",
                CRVAL1 = 0,
                CDELT1 = self.arcsecpix,
                CRPIX1 = self.refpixs[ck][0] + 1,
                CUNIT1 = "arcsec",
                CTYPE2 = "OFFSET",
                CRVAL2 = 0,
                CDELT2 = self.arcsecpix,
                CRPIX2 = self.refpixs[ck][1] + 1,
                CUNIT2 = "arcsec"
            )
            hdu = fits.PrimaryHDU(mom2)
            hdul = fits.HDUList([hdu])
            hdu.header = hdr
            ut.make_folder(foldername=f'models/{self.modelname}/fits')
            hdul.writeto(f'models/{self.modelname}/fits/{ck}_mom2.fits', overwrite=True)
            if self.verbose:
                print(f'models/{self.modelname}/fits/{ck}_mom2.fits saved')
        return mom2

    def mom8(self, ck, chan_range=None, clipping=0, save=False):
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        # chan_vels = self.velchans[chan_range[0]:chan_range[-1]]
        cube_clipped = np.copy(self.cubes[ck])
        clipping = clipping if clipping != 0 \
            else self.momtol_clipping * np.max(self.cubes[ck])
        cube_clipped[cube_clipped<clipping] = 0
        mom8 = np.nan_to_num(
                moments.mom8(
                    cube_clipped,
                    chan_range=chan_range,
                )
            )
        if save:
            hdr = ut.create_hdr(
                NAXIS1 = np.shape(self.cubes[ck])[1],
                NAXIS2 = np.shape(self.cubes[ck])[0],
                BMAJ = self.bmaj/3600,
                BMIN = self.bmin/3600,
                BPA = self.pabeam,
                BTYPE = self.btypes[self.q(ck)],
                BUNIT = self.bunits[self.q(ck)],
                CTYPE1 = "OFFSET",
                CRVAL1 = 0,
                CDELT1 = self.arcsecpix,
                CRPIX1 = self.refpixs[ck][0] + 1,
                CUNIT1 = "arcsec",
                CTYPE2 = "OFFSET",
                CRVAL2 = 0,
                CDELT2 = self.arcsecpix,
                CRPIX2 = self.refpixs[ck][1] + 1,
                CUNIT2 = "arcsec"
            )
            hdu = fits.PrimaryHDU(mom8)
            hdul = fits.HDUList([hdu])
            hdu.header = hdr
            ut.make_folder(foldername=f'models/{self.modelname}/fits')
            hdul.writeto(f'models/{self.modelname}/fits/{ck}_mom8.fits', overwrite=True)
            if self.verbose:
                print(f'models/{self.modelname}/fits/{ck}_mom8.fits saved')
        return mom8

    def plotpv(self, pvimage, rangex, chan_vels, **kwargs):
            return ut.plotpv(pvimage, rangex, chan_vels, **kwargs)

    def plotsumint(self, sumint, **kwargs):
            return ut.plotsumint(sumint, **kwargs)

    def plotmom0(self, mom0, **kwargs):
                return ut.plotmom0(mom0, **kwargs)

    def plotmom1(self, mom1, **kwargs):
                return ut.plotmom1(mom1, **kwargs)

    def plotmom2(self, mom2, **kwargs):
                return ut.plotmom2(mom2, **kwargs)

    def plotmom8(self, mom8, **kwargs):
                return ut.plotmom8(mom8, **kwargs)

    def momentsandpv_all(self, **kwargs):
        for ck in self.listmompvs:
            self.momentsandpv(ck, **kwargs)

    def momentsandpv(self, ck, savefits=False, saveplot=False,
                      mom1clipping=0, mom2clipping=0, verbose=True,
                      mom0values={v: None for v in ["vmax", "vcenter", "vmin"]},
                      mom1values={v: None for v in ["vmax", "vcenter", "vmin"]},
                      mom2values={v: None for v in ["vmax", "vcenter", "vmin"]},
                      mom8values={v: None for v in ["vmax", "vcenter", "vmin"]},
                      pvvalues={v: None for v in ["vmax", "vcenter", "vmin"]},
                      ):
        if verbose:
            print(
"""
\nComputing moments and the PV-diagram along the jet axis
"""
            )

        ckpv = ck + "R"
        if ckpv not in self.cubes:
            self.rotate(ck, forpv=True)

        fig = plt.figure(figsize=(22,9))
        gs = GridSpec(
            2, 5,
            width_ratios=[0.9]*4 + [1],
            hspace=0.2,
            wspace=0.2,
        )
        axs = {}
        cbaxs = {}
        gss = {}

        for i, ik in enumerate(["mom0", "mom1", "mom2", "mom8", "pv"]):
            gss[ik] = gs[0,i].subgridspec(
                2, 1,
                height_ratios=[0.05, 1],
                width_ratios=[1],
                hspace=0.05,
            )
            axs[ik] = plt.subplot(gss[ik][1,0])
            cbaxs[ik] = plt.subplot(gss[ik][0,0])

        ak = "mom0"
        mom0 = self.mom0(
            ck,
            chan_range=[0, self.nc],
            save=savefits,
            )
        extent = np.array([
            -(-0.5-self.refpixs[ck][0]),
            -(self.nxs-0.5-self.refpixs[ck][0]),
            (-0.5-self.refpixs[ck][1]),
            (self.nys-0.5-self.refpixs[ck][1]),
            ]) * self.arcsecpix
        self.plotmom0(
            mom0,
            extent=extent,
            interpolation=None,
            ax=axs[ak],
            cbax=cbaxs[ak],
            cbarlabel="Integrated " + self.getunitlabel(ck).rstrip("]") + " km/s]",
            **mom0values,
            )

        ak = "mom1"
        clipping = float(mom1clipping.split("x")[0]) * self.sigma_noises[ck] \
            if mom1clipping !=0 else 0
        mom1 = self.mom1(
                    ck,
                    chan_range=[0, self.nc],
                    save=savefits,
                    clipping=clipping,
                )
        extent = np.array([
            -(-0.5-self.refpixs[ck][0]),
            -(self.nxs-0.5-self.refpixs[ck][0]),
            (-0.5-self.refpixs[ck][1]),
            (self.nys-0.5-self.refpixs[ck][1]),
            ]) * self.arcsecpix
        self.plotmom1(
            mom1,
            extent=extent,
            interpolation=None,
            ax=axs[ak],
            cbax=cbaxs[ak],
            cbarlabel="Mean velocity [km/s]",
            **mom1values,
            )

        ak = "mom2"
        clipping = float(mom2clipping.split("x")[0]) * self.sigma_noises[ck]\
            if mom1clipping !=0 else 0
        mom2 = self.mom2(
                    ck,
                    chan_range=[0, self.nc],
                    save=savefits,
                    clipping=clipping,
                )
        extent = np.array([
            -(-0.5-self.refpixs[ck][0]),
            -(self.nxs-0.5-self.refpixs[ck][0]),
            (-0.5-self.refpixs[ck][1]),
            (self.nys-0.5-self.refpixs[ck][1]),
            ]) * self.arcsecpix
        self.plotmom2(
            mom2,
            extent=extent,
            interpolation=None,
            ax=axs[ak],
            cbax=cbaxs[ak],
            cbarlabel="Velocity dispersion [km$^2$/s$^2$]",
            **mom2values,
            )

        ak = "pv"
        pvimage = self.pvalongz(
            ckpv,
            halfwidth=0,
            save=savefits
            )
        rangex = np.array([
            -0.5-self.refpixs[ckpv][0],
            self.nxs-0.5-self.refpixs[ckpv][0]
            ]) * self.arcsecpix
        chan_vels = self.velchans
        self.plotpv(
            pvimage,
            rangex=rangex,
            chan_vels=chan_vels,
            ax=axs[ak],
            cbax=cbaxs[ak],
            cbarlabel=self.getunitlabel(ckpv),
            **pvvalues,
            )

        ak = "mom8"
        mom8 = self.mom8(
                    ck,
                    chan_range=[0, self.nc],
                    save=savefits
                )
        extent = np.array([
            -(-0.5-self.refpixs[ck][0]),
            -(self.nxs-0.5-self.refpixs[ck][0]),
            (-0.5-self.refpixs[ck][1]),
            (self.nys-0.5-self.refpixs[ck][1]),
            ]) * self.arcsecpix
        self.plotmom8(
            mom8,
            extent=extent,
            ax=axs[ak],
            cbax=cbaxs[ak],
            cbarlabel="Peak " + self.getunitlabel(ck),
            **mom8values,
            )
        if saveplot:
            fig.savefig(
                f"models/{self.modelname}/momentsandpv.pdf",
                bbox_inches="tight",
                )

    def momentsandpv_and_params_all(self, bscs,**kwargs):
        for ck in self.listmompvs:
            self.momentsandpv_and_params(ck, bscs, **kwargs)

    def momentsandpv_and_params(self, ck, bscs, savefits=False, saveplot=False,
                      mom1clipping=0, mom2clipping=0, verbose=True,
                      mom0values={v: None for v in ["vmax", "vcenter", "vmin"]},
                      mom1values={v: None for v in ["vmax", "vcenter", "vmin"]},
                      mom2values={v: None for v in ["vmax", "vcenter", "vmin"]},
                      mom8values={v: None for v in ["vmax", "vcenter", "vmin"]},
                      pvvalues={v: None for v in ["vmax", "vcenter", "vmin"]},
                      ):
        if verbose:
            print(
"""
\nComputing moments and the PV-diagram along the jet axis
"""
            )


        ckpv = ck + "R"
        if ckpv not in self.cubes:
            self.rotate(ck, forpv=True)

        fig = plt.figure(figsize=(14,10))
        gs = GridSpec(
            2, 3,
            width_ratios=[1]*3, # + [0.85]*2,
            height_ratios=[1]*2,
            hspace=0.3,
            wspace=0.25,
        )
        axs = {}
        cbaxs = {}
        gss = {}

        ik = "text"
        axs[ik] = plt.subplot(gs[0,0])
        axs[ik].set_axis_off()

        ik = "mom0"
        gss[ik] = gs[0,1].subgridspec(
                 2, 1,
                 height_ratios=[0.05, 1],
                 width_ratios=[1],
                 hspace=0.05,
             )
        axs[ik] = plt.subplot(gss[ik][1,0])
        cbaxs[ik] = plt.subplot(gss[ik][0,0])

        ik = "mom8"
        gss[ik] = gs[0,2].subgridspec(
                 2, 1,
                 height_ratios=[0.05, 1],
                 width_ratios=[1],
                 hspace=0.05,
             )
        axs[ik] = plt.subplot(gss[ik][1,0])
        cbaxs[ik] = plt.subplot(gss[ik][0,0])

        ik = "pv"
        gss[ik] = gs[1,0].subgridspec(
                 2, 1,
                 height_ratios=[0.05, 1],
                 width_ratios=[1],
                 hspace=0.05,
             )
        axs[ik] = plt.subplot(gss[ik][1,0])
        cbaxs[ik] = plt.subplot(gss[ik][0,0])


        ik = "mom1"
        gss[ik] = gs[1,1].subgridspec(
                 2, 1,
                 height_ratios=[0.05, 1],
                 width_ratios=[1],
                 hspace=0.05,
             )
        axs[ik] = plt.subplot(gss[ik][1,0])
        cbaxs[ik] = plt.subplot(gss[ik][0,0])

        ik = "mom2"
        gss[ik] = gs[1,2].subgridspec(
                 2, 1,
                 height_ratios=[0.05, 1],
                 width_ratios=[1],
                 hspace=0.05,
             )
        axs[ik] = plt.subplot(gss[ik][1,0])
        cbaxs[ik] = plt.subplot(gss[ik][0,0])

        ak = "text"
        ies = [bsc.i*180/np.pi for bsc in bscs]
        #rjs = [bsc.rj for bsc in bscs]
        L0s = [bsc.L0_arcsec for bsc in bscs]
        zjs = [bsc.zj_arcsec for bsc in bscs]
        vjs = [bsc.vj for bsc in bscs]
        vas = [bsc.va for bsc in bscs]
        v0s = [bsc.v0 for bsc in bscs]
        rbfs = [bsc.rbf_arcsec for bsc in bscs]
        tjs = [bsc.tj_yr for bsc in bscs]
        masss = [bsc.mass*10**4 for bsc in bscs]
        rhows = [bsc.rhow_gcm3*10**20 for bsc in bscs]
        m0s = [bsc.mp0_solmassyr*10**6 for bsc in bscs]
        mwfs = [bsc.mpamb_f_solmassyr*10**6 for bsc in bscs]

        showtext = \
        fr"""
        {bscs[0].modelname}
        Number of bowshocks: {len(bscs)}
        Tex = {self.Tex.value} K
        $i = {{{ut.list2str(ies)}}}^\circ$
        $v_\mathrm{{sys}} = {self.vsys}$ km/s
        $v_\mathrm{{iws}} = {{{ut.list2str(vjs)}}}$ km/s
        $v_0 = {{{ut.list2str(v0s)}}}$ km/s
        $v_a = {{{ut.list2str(vas)}}}$ km/s
        $L_0 = {{{ut.list2str(L0s)}}}$ arcsec
        $z_\mathrm{{iws}} = {{{ut.list2str(zjs)}}}$ arcsec
        $r_\mathrm{{b,f}} = {{{ut.list2str(rbfs)}}}$ arcsec
        $t_\mathrm{{iws}} = {{{ut.list2str(tjs)}}}$ yr
        mass $= {{{ut.list2str(masss)}}}\times 10^{{-4}}$ M$_\odot$
        $\rho_a = {{{ut.list2str(rhows)}}}\times 10^{{-20}}$ g cm$^{{-3}}$
        $\dot{{m}}_0 = {{{ut.list2str(m0s)}}}\times10^{{-6}}$ M$_\odot$ yr$^{{-1}}$
        $\dot{{m}}_{{a,f}} = {{{ut.list2str(mwfs)}}}\times10^{{-6}}$ M$_\odot$ yr$^{{-1}}$
        """
        for n, line in enumerate(showtext.split("\n")):
            axs["text"].text(0, 0.99-0.06*n, line, fontsize=12-len(bscs),
                              transform=axs["text"].transAxes)

        ak = "mom0"
        mom0 = self.mom0(
            ck,
            chan_range=[0, self.nc],
            save=savefits,
            )
        extent = np.array([
            -(-0.5-self.refpixs[ck][0]),
            -(self.nxs-0.5-self.refpixs[ck][0]),
            (-0.5-self.refpixs[ck][1]),
            (self.nys-0.5-self.refpixs[ck][1]),
            ]) * self.arcsecpix
        self.plotmom0(
            mom0,
            extent=extent,
            interpolation=None,
            ax=axs[ak],
            cbax=cbaxs[ak],
            cbarlabel="Integrated " + self.getunitlabel(ck).rstrip("]") + " km/s]",
            **mom0values,
            )

        ak = "mom1"
        clipping = float(mom1clipping.split("x")[0]) * self.sigma_noises[ck] \
            if mom1clipping !=0 else 0
        mom1 = self.mom1(
                    ck,
                    chan_range=[0, self.nc],
                    save=savefits,
                    clipping=clipping,
                )
        extent = np.array([
            -(-0.5-self.refpixs[ck][0]),
            -(self.nxs-0.5-self.refpixs[ck][0]),
            (-0.5-self.refpixs[ck][1]),
            (self.nys-0.5-self.refpixs[ck][1]),
            ]) * self.arcsecpix
        self.plotmom1(
            mom1,
            extent=extent,
            interpolation=None,
            ax=axs[ak],
            cbax=cbaxs[ak],
            cbarlabel="Mean velocity [km/s]",
            **mom1values,
            )

        ak = "mom2"
        clipping = float(mom2clipping.split("x")[0]) * self.sigma_noises[ck]\
            if mom1clipping !=0 else 0
        mom2 = self.mom2(
                    ck,
                    chan_range=[0, self.nc],
                    save=savefits,
                    clipping=clipping,
                )
        extent = np.array([
            -(-0.5-self.refpixs[ck][0]),
            -(self.nxs-0.5-self.refpixs[ck][0]),
            (-0.5-self.refpixs[ck][1]),
            (self.nys-0.5-self.refpixs[ck][1]),
            ]) * self.arcsecpix
        self.plotmom2(
            mom2,
            extent=extent,
            interpolation=None,
            ax=axs[ak],
            cbax=cbaxs[ak],
            cbarlabel="Velocity dispersion [km$^2$/s$^2$]",
            **mom2values,
            )

        ak = "pv"
        pvimage = self.pvalongz(
            ckpv,
            halfwidth=0,
            save=savefits,
            )
        rangex = np.array([
            -0.5-self.refpixs[ckpv][0],
            self.nxs-0.5-self.refpixs[ckpv][0]
            ]) * self.arcsecpix
        chan_vels = self.velchans
        self.plotpv(
            pvimage,
            rangex=rangex,
            chan_vels=chan_vels,
            ax=axs[ak],
            cbax=cbaxs[ak],
            cbarlabel=self.getunitlabel(ckpv),
            **pvvalues,
            )

        ak = "mom8"
        mom8 = self.mom8(
                    ck,
                    chan_range=[0, self.nc],
                    save=savefits,
                )
        extent = np.array([
            -(-0.5-self.refpixs[ck][0]),
            -(self.nxs-0.5-self.refpixs[ck][0]),
            (-0.5-self.refpixs[ck][1]),
            (self.nys-0.5-self.refpixs[ck][1]),
            ]) * self.arcsecpix
        self.plotmom8(
            mom8,
            extent=extent,
            ax=axs[ak],
            cbax=cbaxs[ak],
            cbarlabel="Peak " + self.getunitlabel(ck),
            **mom8values,
            )
        if saveplot:
            fig.savefig(
                f"models/{self.modelname}/momentsandpv_and_params_{ck}.pdf",
                bbox_inches="tight",
                )

