import numpy as np

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm, ListedColormap
from matplotlib import colormaps



from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve

from pathlib import Path

import os

from datetime import datetime

import subprocess
result = subprocess.run(['whoami'], stdout=subprocess.PIPE)
user = result.stdout.decode().strip()

import bowshockpy.moments as moments

def progressbar_bowshock(
        iteration, total, timelapsed, intervaltime,
        decimals = 1, length = 100, fill = '─', printend = "\r"
        ):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + ')' + ' ' * (length - filledLength)
    print(f'  0{bar}{percent}% | {timelapsed:.0f}/{intervaltime*total:.0f}s', end = printend)
    if iteration == total:
        print()

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

def plotpv(pvimage, rangex, chan_vels, ax=None, cbax=None,
        vmax=None, vcenter=None, vmin=None,
        cmap="nipy_spectral", interpolation="bilinear", cbarlabel="Intensity [Jy/beam]",
        ):
    if ax is None or cbax is None:
        plt.figure(figsize=(5,5))
        gs = GridSpec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.1,
            wspace=0.00,
        )
        ax = plt.subplot(gs[1,0])
        cbax = plt.subplot(gs[0, 0])
    else:
        pass
    vmax = vmax if vmax is not None else np.max(pvimage[~np.isnan(pvimage)])
    vmin = vmin if vmin is not None else np.min(pvimage[~np.isnan(pvimage)])
    vcenter = vcenter if vcenter is not None else (vmax - vmin) / 2 + vmin
    norm = TwoSlopeNorm(vmax=vmax, vcenter=vcenter, vmin=vmin)
    chanwidth = chan_vels[1] - chan_vels[0]
    im = ax.imshow(
        pvimage,
        origin="lower",
        extent=[
            rangex[0], rangex[1],
            chan_vels[0]-chanwidth/2,
            chan_vels[-1]+chanwidth/2
            ],
        norm=norm,
        cmap=cmap,
        interpolation=interpolation,
        )
    ax.set_aspect(np.abs(rangex[0]-rangex[-1]) / np.abs(chan_vels[0]-chan_vels[-1]) )
    ax.set_ylabel("Velocity [km/s]")
    ax.set_xlabel("Distance [arcsec]")
    ax.minorticks_on()
    ax.tick_params(
        which="both",
        direction="in",
        top=True,
        right=True,
        color="w",
    )
    plt.colorbar(im, cax=cbax, orientation="horizontal", label=cbarlabel)
    cbax.tick_params(
        axis="x", top=True, bottom=False,
        labelbottom=False, labeltop=True,
        direction="in",
        )
    cbax.set_xlabel(cbarlabel)
    cbax.xaxis.set_label_position("top")

def plotsumint(sumint, ax=None, cbax=None, extent=None,
               vmax=None, vcenter=None, vmin=None,
               interpolation="bilinear", cbarlabel="Intensity", ):

    if ax is None or cbax is None:
        plt.figure(figsize=(5,5.5))
        gs = GridSpec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.1,
            wspace=0.00,
        )
        ax = plt.subplot(gs[1, 0])
        cbax = plt.subplot(gs[0, 0])
    else:
        pass
    vmax = vmax if vmax is not None else np.max(sumint[~np.isnan(sumint)])
    vmin = vmin if vmin is not None else np.min(sumint[~np.isnan(sumint)])
    vcenter = vcenter if vcenter is not None else (vmax - vmin) / 2 + vmin
    norm = TwoSlopeNorm(vmax=vmax, vcenter=vcenter, vmin=vmin)
    im = ax.imshow(
        sumint,
        origin="lower",
        cmap="inferno",
        norm=norm,
        interpolation=interpolation,
        extent=extent,
        )
    if extent is None:
        ax.set_ylabel("Dec. [pixel]")
        ax.set_xlabel("R.A. [pixel]")
    else:
        ax.set_ylabel("Dec. [arcsec]")
        ax.set_xlabel("R.A. [arcsec]")
    ax.minorticks_on()
    ax.tick_params(
        which="both",
        direction="in",
        top=True,
        right=True,
        color="w",
    )
    ax.set_aspect("equal")
    plt.colorbar(im, cax=cbax, orientation="horizontal")
    cbax.tick_params(
        axis="x", top=True, bottom=False,
        labelbottom=False, labeltop=True,
        direction="in",
        )
    cbax.set_xlabel(rf"$\sum\mathrm{{{cbarlabel}}}_i$")
    cbax.xaxis.set_label_position("top")

def plotmom0(mom0, ax=None, cbax=None, extent=None,
            vmax=None, vcenter=None, vmin=None,
            interpolation="bilinear", cbarlabel="Moment 0 [Jy/beam km/s]",):

    if ax is None or cbax is None:
        plt.figure(figsize=(5,5.5))
        gs = GridSpec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.1,
            wspace=0.00,
        )
        ax = plt.subplot(gs[1, 0])
        cbax = plt.subplot(gs[0, 0])
    else:
        pass

    vmax = vmax if vmax is not None else np.max(mom0[~np.isnan(mom0)])
    vmin = vmin if vmin is not None else np.min(mom0[~np.isnan(mom0)])
    vcenter = vcenter if vcenter is not None else (vmax - vmin) / 2 + vmin
    norm = TwoSlopeNorm(vmax=vmax, vcenter=vcenter, vmin=vmin)
    im = ax.imshow(
        mom0,
        origin="lower",
        cmap="inferno",
        norm=norm,
        interpolation=interpolation,
        extent=extent,
        )
    if extent is None:
        ax.set_ylabel("Dec. [pixel]")
        ax.set_xlabel("R.A. [pixel]")
    else:
        ax.set_ylabel("Dec. [arcsec]")
        ax.set_xlabel("R.A. [arcsec]")
    ax.minorticks_on()
    ax.tick_params(
        which="both",
        direction="in",
        top=True,
        right=True,
        color="w",
    )
    ax.set_aspect("equal")
    plt.colorbar(im, cax=cbax, orientation="horizontal", label=cbarlabel)
    cbax.tick_params(
        axis="x", top=True, bottom=False,
        labelbottom=False, labeltop=True,
        direction="in",
        )
    cbax.xaxis.set_label_position("top")


def plotmom1(mom1, ax=None, cbax=None, extent=None,
              vmin=None, vmax=None, vcenter=None,
              extend_cbar="max", return_velcmap=False,
              bg="black", cmap_ref='jet_r',
              interpolation="bilinear", cbarlabel="Moment 1 [km/s]"):
    """
    Moment 1
    """
    if ax is None or cbax is None:
        plt.figure(figsize=(5,5.5))
        gs = GridSpec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.1,
            wspace=0.00,
        )
        ax = plt.subplot(gs[1, 0])
        cbax = plt.subplot(gs[0, 0])
    else:
        pass

    if type(cmap_ref) is str:
        cmap = colormaps[cmap_ref]
    else:
        cmap = cmap_ref
    velcolors = cmap(np.linspace(0, 1, 256))
    if bg == "black":
        bgcolor = np.array([0/256, 0/256, 0/256, 1])
    elif bg == "white":
        bgcolor = np.array([256/256, 256/256, 256/256, 1])
    velcolors[:1, :] = bgcolor

    velcmap = ListedColormap(velcolors)

    if extend_cbar == "max":
        velcmap = ListedColormap(velcolors[::-1])

    vmin = vmin if vmin is not None else np.min(mom1[~np.isnan(mom1)])
    vmax = vmin if vmax is not None else np.max(mom1[~np.isnan(mom1)])
    vcenter = vcenter if vcenter is not None else (vmax - vmin) / 2 + vmin
    norm = TwoSlopeNorm(vcenter=vcenter, vmax=vmax, vmin=vmin)
    im = ax.imshow(
        mom1,
        origin="lower",
        extent=extent,
        norm=norm,
        cmap=velcmap,
        interpolation=interpolation,
    )
    ax.minorticks_on()
    ax.tick_params(
        which="both",
        direction="in",
        top=True,
        right=True,
        color="w",
    )
    if extent is None:
        ax.set_ylabel("Dec. [pixel]")
        ax.set_xlabel("R.A. [pixel]")
    else:
        ax.set_ylabel("Dec. [arcsec]")
        ax.set_xlabel("R.A. [arcsec]")
    ax.set_aspect("equal")
    plt.colorbar(im, cax=cbax, orientation="horizontal",
                 extend=extend_cbar, label=cbarlabel)
    cbax.tick_params(axis="x", top=True, bottom=False, labelbottom=False,
                     labeltop=True, direction="in")
    cbax.xaxis.set_label_position("top")

    if return_velcmap:
        return velcmap
    else:
        pass

def plotmom2(mom2, ax=None, cbax=None, extent=None,
              vmin=None, vmax=None, vcenter=None,
              extend_cbar="max", return_velcmap=False,
              bg="black", cmap_ref='jet_r', cbarlabel="Moment 2 [km$^2$/s$^2$]",
              interpolation=None):
    """
    Moment 1
    """
    if ax is None or cbax is None:
        plt.figure(figsize=(5,5.5))
        gs = GridSpec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.1,
            wspace=0.00,
        )
        ax = plt.subplot(gs[1, 0])
        cbax = plt.subplot(gs[0, 0])
    else:
        pass

    if type(cmap_ref) is str:
        cmap = colormaps[cmap_ref]
    else:
        cmap = cmap_ref
    velcolors = cmap(np.linspace(0, 1, 256))
    if bg == "black":
        bgcolor = np.array([0/256, 0/256, 0/256, 1])
    elif bg == "white":
        bgcolor = np.array([256/256, 256/256, 256/256, 1])
    velcolors[:1, :] = bgcolor
    velcmap = ListedColormap(velcolors)
    if extend_cbar == "max":
        velcmap = ListedColormap(velcolors[::-1])
    vmin = vmin if vmin is not None else np.min(mom2[~np.isnan(mom2)])
    vmax = vmin if vmax is not None else np.max(mom2[~np.isnan(mom2)])
    vcenter = vcenter if vcenter is not None else (vmax - vmin) / 2 + vmin
    norm = TwoSlopeNorm(vcenter=vcenter, vmax=vmax, vmin=vmin)
    im = ax.imshow(
        mom2,
        origin="lower",
        extent=extent,
        norm=norm,
        cmap=velcmap,
        interpolation=interpolation,
    )
    ax.minorticks_on()
    ax.tick_params(
        which="both",
        direction="in",
        top=True,
        right=True,
        color="w",
    )
    if extent is None:
        ax.set_ylabel("Dec. [pixel]")
        ax.set_xlabel("R.A. [pixel]")
    else:
        ax.set_ylabel("Dec. [arcsec]")
        ax.set_xlabel("R.A. [arcsec]")
    ax.set_aspect("equal")
    plt.colorbar(im, cax=cbax, orientation="horizontal",
                extend=extend_cbar, label=cbarlabel)
    cbax.tick_params(axis="x", top=True, bottom=False, labelbottom=False,
                     labeltop=True, direction="in")
    cbax.xaxis.set_label_position("top")
    if return_velcmap:
        return velcmap
    else:
        pass

def plotmom8(mom8, ax=None, cbax=None, extent=None,
            vmax=None, vcenter=None, vmin=None,
            interpolation="bilinear", cbarlabel="Moment 8"):

    if ax is None or cbax is None:
        plt.figure(figsize=(5,5.5))
        gs = GridSpec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.1,
            wspace=0.00,
        )
        ax = plt.subplot(gs[1, 0])
        cbax = plt.subplot(gs[0, 0])
    else:
        pass

    vmax = vmax if vmax is not None else np.max(mom8[~np.isnan(mom8)])
    vmin = vmin if vmin is not None else np.min(mom8[~np.isnan(mom8)])
    vcenter = vcenter if vcenter is not None else (vmax - vmin) / 2 + vmin
    norm = TwoSlopeNorm(vmax=vmax, vcenter=vcenter, vmin=vmin)
    im = ax.imshow(
        mom8,
        origin="lower",
        cmap="inferno",
        norm=norm,
        interpolation=interpolation,
        extent=extent,
        )
    if extent is None:
        ax.set_ylabel("Dec. [pixel]")
        ax.set_xlabel("R.A. [pixel]")
    else:
        ax.set_ylabel("Dec. [arcsec]")
        ax.set_xlabel("R.A. [arcsec]")
    ax.minorticks_on()
    ax.tick_params(
        which="both",
        direction="in",
        top=True,
        right=True,
        color="w",
    )
    ax.set_aspect("equal")
    plt.colorbar(im, cax=cbax, orientation="horizontal", label=cbarlabel)
    cbax.tick_params(
        axis="x", top=True, bottom=False,
        labelbottom=False, labeltop=True,
        direction="in",
        )
    cbax.xaxis.set_label_position("top")

