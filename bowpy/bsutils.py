"""
SOME THINGS FROM THIS FILE IS INTENDED TO DISAPPEAR
Usefull things from bs.py should be moved here
"""

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, ListedColormap
from matplotlib import cm

from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve

from pathlib import Path

import os

from datetime import datetime

import subprocess
result = subprocess.run(['whoami'], stdout=subprocess.PIPE)
user = result.stdout.decode().strip()

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