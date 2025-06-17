from matplotlib import colormaps
from matplotlib import colors

import numpy as np

from astropy import units as u

import subprocess

process = subprocess.Popen(["whoami"], stdout=subprocess.PIPE)
result = process.communicate()[0]
user = result.decode('utf-8').rstrip('\n')

def list2str(a, precision=2):
    _list = [float(f'{i:.{precision}f}') for i in a]
    _str = str(_list) if len(_list)>1 else str(_list[0])
    return _str

def get_color(vel_range, vel, cmap, norm="linear"):
    """
    Gets the color that corresponds in a colormap linearly interpolated taking
    into account the values at the limits.
    """
    cmapp = colormaps.get_cmap(cmap)
    if norm == "linear":
        norm = colors.Normalize(vmin=vel_range[0], vmax=vel_range[-1])
    elif norm == "log":
        norm = colors.LogNorm(vmin=vel_range[0], vmax=vel_range[-1])
    rgba = cmapp(norm(vel))
    color = colors.to_hex(rgba)
    return color

def mb_sa_gaussian_f(maja, mina):
    """
    Solid angle of a gaussian main beam and θmaj and θmin as
    the half-power beam widths
    """
    omega_M = np.pi * maja * mina / (4 * np.log(2))
    return omega_M.to(u.sr)
