import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, Normalize, ListedColormap
from matplotlib.widgets import Slider, Button, TextBox
from matplotlib.gridspec import GridSpec
from matplotlib import cm

from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.coordinates import SkyCoord, concatenate
import astropy.units as u
from astropy.wcs import WCS

from scipy.optimize import minimize, minimize_scalar
from scipy.integrate import quad

import shelve

from pathlib import Path

import os

from datetime import datetime

import subprocess
result = subprocess.run(['whoami'], stdout=subprocess.PIPE)
user = result.stdout.decode().strip()

import bowpy.utils as utils
import bowpy.moments as moments


path = Path(__file__).parent / "header_default.txt"
with path.open() as f:
    hdr_str_default = f.read()

def create_hdr(**kwargs):
    hdr = fits.Header.fromstring(hdr_str_default)
    if len(kwargs) != 0:
         for kwarg in kwargs:
             hdr[kwarg] = kwargs[kwarg]
    return hdr

def make_folder(foldername=None):
    if not os.path.exists(foldername):
        os.makedirs(foldername)

def write_log(path, mode, pars, printtime=False):
    with open(path, mode) as f:
        if printtime:
            time = datetime.now().strftime("%d/%m/%y %H:%M:%S")
            f.writelines(f"###################\n# {time}\n###################\n")
        f.writelines("pars = {\n")
        for par in pars:
            f.writelines(f'    "{par}": "{pars[par]}",\n')
        f.writelines("}\n\n")

def read_modelparams(modelname, namepars="params"):
    """
    Reads the params.txt within the folder modelname
    """
    pars = {}
    path_bs = "models"
    with open(f"{path_bs}/{modelname}/{namepars}.txt") as f:
        for i,line in enumerate(f):
            if ("{" in line) or ("}" in line) or (line=="\n"):
                pass
            else:
                ppp = line
                pars[line.split(':')[0].strip().strip('"')] = \
                        line.split(":")[1].rstrip(",\n").strip().strip('"')
    for par in pars:
        try:
            pars[par] = float(pars[par])
        except:
            pass
    return pars

def gaussconvolve(data, x_FWHM, y_FWHM, pa, return_kernel=False):
    """
    Gausskernel 0 and 1 entries are the FWHM, the third the PA
    """
    x_stddev = x_FWHM / (2 * np.sqrt(2 * np.log(2)))
    y_stddev = y_FWHM / (2 * np.sqrt(2 * np.log(2)))
    kernel = Gaussian2DKernel(
        x_stddev=x_stddev,
        y_stddev=y_stddev,
        theta=pa*np.pi/180)
    data_conv = convolve(data, kernel)
    if return_kernel:
        return data_conv, kernel
    else:
        return data_conv

def channels_plot(cube, axs, ax_cbar, pars, chan_vels, nrow, ncol,
                  fmaxlim=0.15, fvcenter=0.1, vmax=None, vcenter=None,
                  vmin=None, plotparams=False):
    """
    Plots some channels from bowshock_cube. The number of channels is
    determined by nrow*ncol
    """
    chans_plot = float(pars["NC"])
    nv = np.shape(cube)[2]
    vels = chan_vels

    if plotparams:
        nchanscube = nrow*ncol-1
    else:
        nchanscube = nrow*ncol
    selint = int(chans_plot/nchanscube)
    alldata = cube[::selint]
    if vmax is None:
        uplim = np.max(alldata) * fmaxlim
        norm = TwoSlopeNorm(vcenter=uplim*fvcenter, vmax=uplim, vmin=0)
    else:
        vmin = vmin if vmin is not None else 0
        vcenter = vcenter if vcenter is not None else (vmax - vmin) / 2.
        norm = TwoSlopeNorm(vmax=vmax, vcenter=vcenter, vmin=vmin)

    for cha in range(nchanscube):
        if plotparams:
            chan = cha + 1
        else:
            chan = cha
        data = cube[::selint][chan]
        im = axs[chan].imshow(
            data,
            norm=norm,
            origin="lower", interpolation="bilinear", cmap="inferno",)

        axs[chan].text(0.05, 0.9,
                s=f"V={chan_vels[::selint][chan]:.2f}",
                color="w", transform=axs[chan].transAxes, fontsize=15)

        axs[chan].set_xlabel("Ra")
        axs[chan].set_ylabel("Dec")
        axs[chan].set_aspect("equal")
    plt.colorbar(im, cax=ax_cbar, extend="max")
    ax_cbar.set_ylabel("Intensity")

def pv_plot(cube, ax, cbax, chan_vels, rangex=None, width_pv=3, xpv=None,
            fmaxlim=0.4, fvcenter=0.2, vmax=None, vcenter=None, vmin=None,
            cmap="inferno", interpolation="bilinear", normalize=False):
    """
    PV
    """
    xpv = xpv if xpv is not None else int(np.shape(cube)[-1]/2)
    bowshock_pv = moments.pv(cube, xpv=xpv, width=width_pv)
    normfactor = np.max(bowshock_pv[:, ::-1]) if normalize else 1
    data = bowshock_pv[:, ::-1] / normfactor
    rangex = rangex if rangex is not None else [-0.5, np.shape(data)[1]-0.5]
    if vmax is None:
        maxlim = np.max(data) * fmaxlim
        norm = TwoSlopeNorm(vcenter=maxlim*fvcenter, vmax=maxlim, vmin=0)
    else:
        vmin = vmin if vmin is not None else 0
        vcenter = vcenter if vcenter is not None else (vmax - vmin) / 2.
        norm = TwoSlopeNorm(vmax=vmax, vcenter=vcenter, vmin=vmin)


    im = ax.imshow(
        data,
        origin="lower",
        extent=[rangex[0], rangex[1],
                chan_vels[0]-np.abs(chan_vels[0]-chan_vels[1])/2.,
                chan_vels[-1]-np.abs(chan_vels[0]-chan_vels[1])/2.],
        norm=norm,
        cmap=cmap,
        interpolation=interpolation)

    plt.colorbar(im, cax=cbax, orientation="horizontal", extend="max")
    cbax.tick_params(axis="x", top=True, bottom=False, labelbottom=False,
                     labeltop=True,)
    cbax.set_xlabel("Intensity")
    cbax.xaxis.set_label_position("top")
    ax.set_aspect("auto")
    ax.set_ylabel("Velocity (km/s)")
    ax.set_xlabel("Distance (arcsec)")
    return data

def sumint_plot(cube, ax, cbax, pars, chan0=None, chanf=None,
                fmaxlim=0.4, fvcenter=0.2, vmax=None, vcenter=None,
                vmin=None):
    """
    Sumint
    """
    chan0 = chan0 if chan0 is not None else 0
    chanf = chanf if chanf is not None else np.shape(cube)[0]
    data = moments.sumint(cube,
                         chan_range=[chan0, chanf])
    if vmax is None:
        uplim = np.max(data) * fmaxlim
        norm = TwoSlopeNorm(vcenter=uplim*fvcenter, vmax=uplim, vmin=0)
    else:
        vmin = vmin if vmin is not None else 0
        vcenter = vcenter if vcenter is not None else (vmax - vmin) / 2.
        norm = TwoSlopeNorm(vmax=vmax, vcenter=vcenter, vmin=vmin)
    im = ax.imshow(data,
              origin="lower",
              cmap="inferno",
              norm=norm,
              interpolation="bilinear")

    plt.colorbar(im, cax=cbax, orientation="horizontal", extend="max",)
    cbax.tick_params(axis="x", top=True, bottom=False, labelbottom=False,
                     labeltop=True,)
    cbax.set_xlabel(r"$\sum\mathrm{I}_i$", labelpad=10, fontsize=12)
    cbax.xaxis.set_label_position("top")
    ax.set_aspect("auto")
    ax.set_ylabel("Dec (pixel)")
    ax.set_xlabel("Ra (pixel)")

    ax.set_aspect("equal")
    return data

def mom0_plot(cube, ax, cbax, pars, chan_vels,
              chan0=None, chanf=None, fmaxlim=0.4, fvcenter=0.2,
              vmax=None, vcenter=None, vmin=None):
    """
    Moment 0
    """
    chan0 = chan0 if chan0 is not None else 0
    chanf = chanf if chanf is not None else np.shape(cube)[0]

    data = moments.mom0(cube,
                     chan_vels=chan_vels,
                     chan_range=[chan0, chanf])
    if vmax is None:
        uplim = np.max(data) * fmaxlim
        norm = TwoSlopeNorm(vcenter=uplim*fvcenter, vmax=uplim, vmin=0)
    else:
        vmin = vmin if vmin is not None else 0
        vcenter = vcenter if vcenter is not None else (vmax - vmin) / 2.
        norm = TwoSlopeNorm(vmax=vmax, vcenter=vcenter, vmin=vmin)
    im = ax.imshow(data,
              origin="lower",
              cmap="inferno",
              norm=norm,
              interpolation="bilinear")
    plt.colorbar(im, cax=cbax, orientation="horizontal", extend="max",)
    cbax.tick_params(axis="x", top=True, bottom=False, labelbottom=False,
                     labeltop=True,)
    cbax.set_xlabel(r"Moment 0", labelpad=10, fontsize=12)
    cbax.xaxis.set_label_position("top")
    ax.set_aspect("auto")
    ax.set_ylabel("Dec (pixel)")
    ax.set_xlabel("Ra (pixel)")

    ax.set_aspect("equal")
    return data

def mom1_plot(cube, ax, cbax, pars, chan_vels,
              chan0=None, chanf=None, vmin=-100, vmax=20, vcenter=-50,
              extend_cbar="max", clipping=0, return_velcmap=False,
              show_plot=True, bg="black", cmap_ref='jet_r',
              interpolation=None):
    """
    Moment 1
    """
    if type(cmap_ref) is str:
        cmap = cm.get_cmap(cmap_ref, 256)
    else:
        cmap = cmap_ref
    velcolors = cmap(np.linspace(0, 1, 256))
    if bg == "black":
        bgcolor = np.array([0/256, 0/256, 0/256, 1])
    elif bg == "white":
        bgcolor = np.array([256/256, 256/256, 256/256, 1])
    velcolors[:1, :] = bgcolor

    velcmap = ListedColormap(velcolors)

    chan0 = chan0 if chan0 is not None else 0
    chanf = chanf if chanf is not None else np.shape(cube)[0]
    cube_clipped = np.copy(cube)
    cube_clipped[cube_clipped<clipping] = 0
    data = np.nan_to_num(
              moments.mom1(
                 cube_clipped,
                 chan_vels=chan_vels,
                 chan_range=[chan0, chanf])
                 )

    if extend_cbar == "max":
        velcmap = ListedColormap(velcolors[::-1])
    if show_plot:
        im = ax.imshow(
            data,
            origin="lower",
            norm=TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax),
            cmap=velcmap,
            interpolation=interpolation,
        )

        plt.colorbar(im, cax=cbax, orientation="horizontal", extend=extend_cbar,)
        cbax.tick_params(axis="x", top=True, bottom=False, labelbottom=False,
                         labeltop=True,)
        cbax.set_xlabel(r"Moment 1 (km/s)", labelpad=10, fontsize=12)
        cbax.xaxis.set_label_position("top")

        ax.set_aspect("auto")
        ax.set_ylabel("Dec (pixel)")
        ax.set_xlabel("Ra (pixel)")

        ax.set_aspect("equal")
    else:
        pass
    if return_velcmap:
        return data, velcmap
    else:
        return data

def mom2_plot(cube, ax, cbax, pars, chan_vels,
              chan0=None, chanf=None, vmin=1, vcenter=10,
              vmax=15, extend_cbar="both", clipping=0):
    """
    Moment 2
    """
    cmap = cm.get_cmap('jet_r', 256)
    velcolors = cmap(np.linspace(0, 1, 256))
    blackcolor = np.array([0/256, 0/256, 0/256, 1])
    velcolors[:1, :] = blackcolor
    velcmap = ListedColormap(velcolors)

    chan0 = chan0 if chan0 is not None else 0
    chanf = chanf if chanf is not None else np.shape(cube)[0]
    cube_clipped = np.copy(cube)
    cube_clipped[cube_clipped<clipping] = 0
    data =  np.nan_to_num(
                moments.mom2(
                    cube_clipped,
                    chan_vels=chan_vels,
                    chan_range=[chan0, chanf])
                    )
    # if data[np.where(np.abs(data)==np.max(np.abs(data)))] >= 0:
    #     uplim = np.max(data)
    #     lowlim = uplim * fminlim
    #     extend_cbar="min"
    # else:
    #     lowlim = np.min(data)
    #     uplim = lowlim * fminlim
    #     velcmap = ListedColormap(velcolors[::-1])
    #     extend_cbar="max"
    im = ax.imshow(
        data,
        norm=TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax),
        origin="lower",
        cmap=velcmap
    )
    plt.colorbar(im, cax=cbax, orientation="horizontal", extend=extend_cbar,)
    cbax.tick_params(axis="x", top=True, bottom=False, labelbottom=False,
                     labeltop=True,)
    cbax.set_xlabel(r"Moment 2 (km/s)", labelpad=10, fontsize=12)
    cbax.xaxis.set_label_position("top")
    ax.set_aspect("auto")
    ax.set_ylabel("Dec (pixel)")
    ax.set_xlabel("Ra (pixel)")

    ax.set_aspect("equal")
    return data

def mom8_plot(cube, ax, cbax, pars, chan0=None, chanf=None,
              fmaxlim=0.4, fvcenter=0.2, vmax=None, vcenter=None,
              vmin=None):
    """
    Moment 8
    """
    chan0 = chan0 if chan0 is not None else 0
    chanf = chanf if chanf is not None else np.shape(cube)[0]
    data = moments.mom8(cube,
                         chan_range=[chan0, chanf])
    if vmax is None:
        uplim = np.max(data) * fmaxlim
        norm = TwoSlopeNorm(vcenter=uplim*fvcenter, vmax=uplim, vmin=0)
    else:
        vmin = vmin if vmin is not None else 0
        vcenter = vcenter if vcenter is not None else (vmax - vmin) / 2.
        norm = TwoSlopeNorm(vmax=vmax, vcenter=vcenter, vmin=vmin)
    im = ax.imshow(data,
              origin="lower",
              norm=norm,
              cmap="inferno")
    plt.colorbar(im, cax=cbax, orientation="horizontal", extend="max",)
    cbax.tick_params(axis="x", top=True, bottom=False, labelbottom=False,
                     labeltop=True,)
    cbax.set_xlabel(r"Moment 8", labelpad=10, fontsize=12)
    cbax.xaxis.set_label_position("top")
    ax.set_aspect("auto")
    ax.set_ylabel("Dec (pixel)")
    ax.set_xlabel("Ra (pixel)")

    ax.set_aspect("equal")
    return data


def moments_sheet(cube, axs, cbaxs, pars, chan_vels,
                  fmaxlim_pv=0.4, fvcenter_pv=0.2,
                  fmaxlim_sumint=0.4, fvcenter_sumint=0.2,
                  fmaxlim_mom0=0.4, fvcenter_mom0=0.2,
                  fminlim_mom1=0.2, fvcenter_mom1=0.6,
                  fminlim_mom2=0.2, fvcenter_mom2=3.2,
                  fmaxlim_mom8=0.4, fvcenter_mom8=0.2):
    """
    Generate a sheet with a PV, mom0, mom1, mom2 and mom8
    """

    """
    PV
    """
    _ = pv_plot(cube, axs[0], cbaxs[0], chan_vels,
                fmaxlim=fmaxlim_pv, fvcenter=fvcenter_pv)

    """
    Intintens
    """
    _ = sumint_plot(cube, axs[1], cbaxs[1], pars,
                    fmaxlim=fmaxlim_sumint, fvcenter=fvcenter_sumint)

    """
    Moment 0
    """
    _ = mom0_plot(cube, axs[2], cbaxs[2], pars, chan_vels,
                    fmaxlim=fmaxlim_mom0, fvcenter=fvcenter_mom0)

    """
    Moment 1
    """
    _ = mom1_plot(cube, axs[3], cbaxs[3], pars, chan_vels,
                    fminlim=fminlim_mom1, fvcenter=fvcenter_mom1)

    """
    Moment 2
    """
    _ = mom2_plot(cube, axs[4], cbaxs[4], pars, chan_vels,
                    fminlim=fminlim_mom2, fvcenter=fvcenter_mom2)

    """
    Moment 8
    """
    _ = mom8_plot(cube, axs[5], cbaxs[5], pars,
                    fmaxlim=fmaxlim_mom8, fvcenter=fvcenter_mom8)


class NJ():

    default_kwargs = {
        "rbf_niters": 1000,
    }

    def __init__(self, ps, **kwargs):
        for param in ps:
            setattr(self, param.replace("/", "_"), ps[param])
        for kwarg in self.default_kwargs:
            kwarg_attr = kwargs[kwarg] if kwarg in kwargs else self.default_kwargs[kwarg]
            setattr(self, kwarg, kwarg_attr)

        self.zj = self.zj if self.zj is not None \
            else (self.xhead - self.xorigin) / np.sin(self.i)
        self.vj = self.vj if self.vj is not None \
            else (-self.vhead + self.vsys) / np.cos(self.i)
        self.tj = (self.zj * mf.default_params["SVS13_distance"] * u.au \
                   / (self.vj * u.km/u.s)).to(u.yr).value
        self.rbf = self.rbf_calc()
        self.zbf = self.zb_r(self.rbf)

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

    def vz_obs(self, zb, phi):
        a = self.alpha(zb)
        return self.vtot(zb) * (np.cos(a)*np.cos(self.i) - np.sin(a)*np.cos(phi)*np.sin(self.i))

    def x_obs(self, zb, phi):
        return self.rb(zb)*np.cos(phi)*np.cos(self.i) + zb*np.sin(self.i)

    def y_obs(self, zb, phi):
        return self.rb(zb) * np.sin(phi)

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
        if self.rhow_u == "g/cm3":
            arcsec2cm = (self.distpc * u.au).to(u.cm).value
            sd = 0.5 * self.rhow * cosa * (self.gamma()*tana+1)**2 * self.rb(zb) * arcsec2cm
        if self.rhow_u == "msun/arcsec3":
            sd = 0.5 * self.rhow * cosa * (self.gamma()*tana+1)**2 * self.rb(zb)
        return sd

    def dr_func(self, zb, dz):
        """
        Differential of r given a differential of z
        """
        return 1/3 * (self.L0 / self.rb(zb))**2 * dz

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
        u = rbf / self.L0 * (3/self.gamma())**0.5
        analit_int = u**5 / 5 + 2*u**3/3 + u
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
        u = rb / self.L0 * (3/self.gamma())**0.5
        analit_int = u**5 / 5 + 2*u**3/3 + u
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
        self.tj = (self.zj * mf.default_params["SVS13_distance"] * u.au \
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


class BowshockFitter():

    limit_slider = {
        'L0_u': 3,
        'L0_l': 0,
        'vw_u': 40,
        'vw_l': 0,
        'v0_u': 60,
        'v0_l': 0,
        'i_u': np.pi/2,
        'i_l': 0,
        'vhead_l': -130,
        'vhead_u': -20,
        'xhead_u': 6,
        'xhead_l': 0,
        'xorigin_u': 7,
        'xorigin_l': -5,
        'rj_u': 1,
        'rj_l': 0,
        }

    nzs_init = 200
    nvisos_init = 200

    def __init__(self, ps, rk, model="nj", **kwargs):
        """
        model: "nj", for narrow jet limit, "aj" for arbitrary width jet
        """
        self.nzs = self.nzs_init if "nzs" not in kwargs else kwargs["nzs"]
        self.nvisos = self.nvisos_init if "nvisos" not in kwargs else kwargs["nvisos"]
        self.kwargs_model = {} if "kwargs_model" not in kwargs \
           else kwargs["kwargs_model"]

        self.param_str = ["L0", "vw", "v0", "i", "vhead", "xhead", "xorigin"]

        self.fig = plt.figure(figsize=(8,8))
        self.fig.subplots_adjust(left=0.25, bottom=0.35)

        self.ps = ps
        if model=="nj":
            self.modclass = NJ
        elif model=="aj":
            self.modclass = AJ
            self.param_str = self.param_str + ["rj"]

        self.mod = self.modclass(self.ps, **self.kwargs_model)

        self.params = {param: getattr(self.mod, param)
                       for param in self.param_str}

        self.rk = rk

        self.fig_model = None
        self.axs = {}
        self.slider_ax = None
        self.create_axes()

        self.param_sliders = None
        self.update_buttons()

        self.plot_model()

    def create_axes(self):
        self.slider_ax = {param: self.fig.add_axes([0.25,0.25-i*0.03,0.65,0.03])
                          for i,param in enumerate(self.params)}

        nrow = 2
        ncol = 6
        wspace = 0.0
        hspace = 0.3
        width_ratios = [0.6, 0.2, 1, 0.02, 0.4, 0.65]
        height_ratios = [1] * nrow

        self.fig_model = plt.figure(figsize=(13,7))
        gs = GridSpec(nrow, ncol,
                      height_ratios=height_ratios,
                      width_ratios=width_ratios)

        gs.update(left=0.05, right=0.95,
                  bottom=0.08, top=0.93,
                  wspace=wspace, hspace=hspace)

        self.axs["text"] = plt.subplot(gs[0, 0])
        self.axs["alpha"] = plt.subplot(gs[1, 0])
        self.axs["cbar"] = plt.subplot(gs[1, 3])
        self.axs[0] = plt.subplot(gs[0, 2:4])
        self.axs[1] = plt.subplot(gs[1, 2:3])
        self.axs[2] = plt.subplot(gs[0, 5])
        self.axs[3] = plt.subplot(gs[1, 5])

        self.ax_text = self.fig.add_axes([0.05,0.025,0.1,0.04])
        self.ax_text.set_axis_off()

    def update_buttons(self,):
        """
        Updates the state of the sliders and buttons (for example, after a fit)
        """

        self.param_sliders = {param: Slider(self.slider_ax[param],
                                            param,
                                            self.limit_slider[param+'_l'],
                                            self.limit_slider[param+'_u'],
                                            valinit=self.params[param])
                              for param in self.params}

        for param in self.params:
            self.param_sliders[param].on_changed(self.sliders_on_changed)

    def update_params(self, params):
        print("Recalculating model...")
        self.params = params
        for param in self.param_str:
            self.ps[param] = params[param]
        self.mod = self.modclass(self.ps)
        self.plot_model()
        print("Updated!\n")

    def plot_model(self,):

        zsextended = self.mod.zb_r(
            np.linspace(self.mod.rbf, self.mod.rj, self.nzs)
        )
        dzs = (np.diff(zsextended[1:]) + np.diff(zsextended[:-1])) / 2
        zs = zsextended[1:-1]

        zs_fromsource = zs + self.mod.xorigin / np.sin(self.mod.i)
        Rs = np.array([self.mod.rb(zb) for zb in zs])
        vrs = np.array([self.mod.vr(zb) for zb in zs])
        vzs = np.array([self.mod.vz(zb) for zb in zs])
        vs = np.array([np.sqrt(vrr**2+vzz**2) for vrr, vzz in zip(vrs, vzs)])
        alphas = np.array([self.mod.alpha(zb) for zb in zs])
        thetas = np.array([np.arctan(Rs[i] / zsource)
                           for i, zsource in enumerate(zs_fromsource)])

        xs_obs90 = np.array([self.mod.x_obs(zb,phi=np.pi/2) for zb in zs]) \
                + self.mod.xorigin
        xs_obs0 = np.array([self.mod.x_obs(zb,phi=0) for zb in zs]) \
                + self.mod.xorigin
        xs_obs180 = np.array([self.mod.x_obs(zb,phi=np.pi) for zb in zs]) \
                + self.mod.xorigin
        Rs_obs = Rs

        vzs_obs0 = -np.array([self.mod.vz_obs(zb,phi=0) for zb in zs]) \
                + self.mod.vsys
        vzs_obs90 = -np.array([self.mod.vz_obs(zb,phi=np.pi/2) for zb in zs]) \
                + self.mod.vsys
        vzs_obs180 = -np.array([self.mod.vz_obs(zb,phi=np.pi) for zb in zs]) \
                + self.mod.vsys

        narrows = 20
        idx_arr = int(len(zs)/narrows)
        larrow = 1 / np.max([np.max(vrs), np.max(vzs)])

        zs_arrows = zs[::idx_arr]
        Rs_arrows = Rs[::idx_arr]

        z_arrow_ref = self.mod.zj
        R_arrow_ref = -0.775

        zs_arrows_tip = zs_arrows + larrow * vzs[::idx_arr] # vs[::idx_arr]*np.cos(alphas[::idx_arr])
        Rs_arrows_tip = Rs_arrows + larrow * vrs[::idx_arr] # vs[::idx_arr]*np.sin(alphas[::idx_arr])

        v_arrow_ref = 100
        z_arrow_ref_tip = z_arrow_ref + larrow * v_arrow_ref
        R_arrow_ref_tip = R_arrow_ref + larrow * 0

        maxviso = np.max(vzs_obs180)
        minviso = np.min(vzs_obs0)
        visos = np.linspace(maxviso, minviso, self.nvisos)

        # Calculate Major and minor axis ratio
        qratios = []
        for viso in visos:
            idx_180 = np.argmin(np.abs((vzs_obs180-viso)))
            idx_90 = np.argmin(np.abs((vzs_obs90-viso)))
            idx_0 = np.argmin(np.abs((vzs_obs0-viso)))

            major_axis = xs_obs0[idx_0] - xs_obs180[idx_180]
            minor_axis = 2 * Rs_obs[idx_90]
            qratios.append(minor_axis / major_axis)

        for axstr in self.axs:
            self.axs[axstr].clear()

        showtext = \
            fr"""
            Rings {self.rk}

            $i = {{{self.mod.i*180/np.pi:.2f}}}^\circ$
            $V_\mathrm{{jet}} = {{{self.mod.vj:.2f}}}$ km/s
            $V_0 = {{{self.mod.v0:.2f}}}$ km/s
            $V_w = {{{self.mod.vw:.2f}}}$ km/s
            $L_0 = {{{self.mod.L0:.2f}}}$ arcsec
            $z_\mathrm{{jet}} = {{{self.mod.zj:.2f}}}$ arcsec
            $r_\mathrm{{b,f}} = {{{self.mod.rbf:.2f}}}$ arcsec
            $t_\mathrm{{jet}} = {{{self.mod.tj:.2f}}}$ yr
            """

        self.axs["text"].set_axis_off()
        for n, line in enumerate(showtext.split("\n")):
             self.axs["text"].text(0, 0.99-0.1*n, line, fontsize=12,
                              transform=self.axs["text"].transAxes)

        """
        Alphas
        """
        self.axs["alpha"].plot(zs, alphas * 180 / np.pi, label="Velocity vector angle")
        self.axs["alpha"].plot(zs, thetas * 180 / np.pi, label="Position vector angle")
        self.axs["alpha"].legend()
        self.axs["alpha"].set_xlabel("$z_b$ (arcsec)")
        self.axs["alpha"].set_ylabel("Angle (degrees)")
        self.axs["alpha"].set_ylim([0, 40])

        """
        Deprojected shell
        """
        self.axs[0].plot(zs, Rs, "b")
        self.axs[0].plot(zs, -Rs, "b")

        for i in range(len(zs_arrows)):
            self.axs[0].annotate(
                "",
                xy=(zs_arrows_tip[i], Rs_arrows_tip[i]),
                xytext=(zs_arrows[i], Rs_arrows[i]),
                arrowprops=dict(arrowstyle="->"))
            self.axs[0].annotate(
                "",
                xy=(zs_arrows_tip[i], -Rs_arrows_tip[i]),
                xytext=(zs_arrows[i], -Rs_arrows[i]),
                arrowprops=dict(arrowstyle="->"))
        self.axs[0].annotate(
            "",
            xy=(z_arrow_ref_tip, R_arrow_ref_tip),
            xytext=(z_arrow_ref, R_arrow_ref),
            arrowprops=dict(arrowstyle="->"))
        self.axs[0].text(z_arrow_ref+0.075,
                    R_arrow_ref+0.05,
                    f"{v_arrow_ref:d} km/s")

        self.axs[0].set_aspect("equal")
        self.axs[0].set_xlabel("$z_b$ (arcsec)")
        self.axs[0].set_ylabel("R (arcsec)")
        self.axs[0].set_xlim([np.min(zs), np.max(zs)*1.2])
        self.axs[0].set_ylim([-np.max(Rs)*1.3, np.max(Rs)*1.3])

        """
        Radius vs distance
        """
        cmap_vel = "jet"
        v_range = [np.min([np.min(visos)]),
                   np.max([np.max(visos)])]
        c_ell = lambda vel: utils.get_color(
            v_range,
            vel,
            cmap_vel
        )
        for i,vel in enumerate(vzs_obs90):
            x = xs_obs90[i]
            y = Rs_obs[i]
            self.axs[1].plot(x, y, ".",
                    color=c_ell(vel),
                    markersize=2)
        self.axs[1].plot(xs_obs90, Rs_obs, "ok", zorder=0)

        self.axs[1].set_aspect("equal")
        self.axs[1].set_xlabel("$x_\mathrm{centers}$ (arcsec)")
        self.axs[1].set_ylabel("R (arcsec)")

        plt.colorbar(cm.ScalarMappable(norm=Normalize(vmax=v_range[1], vmin=v_range[0]),
                                         cmap=cmap_vel,), cax=self.axs["cbar"])
        self.axs["cbar"].invert_yaxis()
        self.axs["cbar"].tick_params(direction='in', color='k', pad=2.)
        self.axs["cbar"].set_ylabel("$V_\mathrm{LSR}$ (km/s)")

        """
        PV diagram
        """
        self.axs[2].plot(xs_obs180, vzs_obs180, "-b")
        self.axs[2].plot(xs_obs0, vzs_obs0, "-b")


        self.axs[2].set_xlabel("$x_\mathrm{obs}$ (arcsec)")
        self.axs[2].set_ylabel("$V_\mathrm{LSR}$ (km/s)")
        self.axs[2].invert_yaxis()

        """
        Axis ratio
        """
        self.axs[3].plot(visos, qratios, "b")
        self.axs[3].invert_xaxis()
        self.axs[3].set_xlabel("$V_\mathrm{LSR}$ (km/s)")
        self.axs[3].set_ylabel("Longitudinal axis / Transversal axis")
        self.axs[3].set_ylim([0.5, 1.1])

        self.fig_model.canvas.draw_idle()

    def sliders_on_changed(self, val):
        self.update_params({param: self.param_sliders[param].val
                            for param in self.param_sliders})
        self.fig.canvas.draw_idle()

