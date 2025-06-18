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

_header_default = "SIMPLE  =                    T / conforms to FITS standard                      BITPIX  =                  -64 / array data type                                NAXIS   =                    3 / number of array dimensions                     NAXIS1  =                  650                                                  NAXIS2  =                  850                                                  NAXIS3  =                  130                                                  BMAJ    =   4.539020773437E-05                                                  BMIN    =   2.322353008721E-05                                                  BPA     =  -2.697280883789E+00                                                  BTYPE   = 'Intensity'                                                           OBJECT  = 'MODEL   '                                                            BUNIT   = 'Jy/beam '           /Brightness (pixel) unit                         RADESYS = 'ICRS    '                                                            LONPOLE =   1.800000000000E+02                                                  LATPOLE =   3.126777777778E+01                                                  PC1_1   =   1.000000000000E+00                                                  PC2_1   =   0.000000000000E+00                                                  PC3_1   =   0.000000000000E+00                                                  PC1_2   =   0.000000000000E+00                                                  PC2_2   =   1.000000000000E+00                                                  PC3_2   =   0.000000000000E+00                                                  PC1_3   =   0.000000000000E+00                                                  PC2_3   =   0.000000000000E+00                                                  PC3_3   =   1.000000000000E+00                                                  CTYPE1  = 'RA---SIN'                                                            CRVAL1  =   5.226562499999E+01                                                  CDELT1  =  -3.333333333333E-06                                                  CRPIX1  =                526.0                                                  CUNIT1  = 'deg     '                                                            CTYPE2  = 'DEC--SIN'                                                            CRVAL2  =   3.126777777778E+01                                                  CDELT2  =   3.333333333333E-06                                                  CRPIX2  =                796.0                                                  CUNIT2  = 'deg     '                                                            CTYPE3  = 'FREQ    '                                                            CRVAL3  =       345821840500.0                                                  CDELT3  =   6.103500000000E+05                                                  CRPIX3  =   1.000000000000E+00                                                  CUNIT3  = 'Hz      '                                                            PV2_1   =   0.000000000000E+00                                                  PV2_2   =   0.000000000000E+00                                                  RESTFRQ =   3.457959900000E+11 /Rest Frequency (Hz)                             SPECSYS = 'LSRK    '           /Spectral reference frame                        ALTRVAL =   8.871029325682E+04 /Alternate frequency reference value             ALTRPIX =   1.000000000000E+00 /Alternate frequency reference pixel             VELREF  =                  257 /1 LSR, 2 HEL, 3 OBS, +256 Radio                 COMMENT casacore non-standard usage: 4 LSD, 5 GEO, 6 SOU, 7 GAL                 TELESCOP= 'ALMA    '                                                            OBSERVER= 'galileo '                                                            DATE-OBS= '2016-09-09T08:27:20.928001'                                          TIMESYS = 'UTC     '                                                            OBSRA   =   5.226562499999E+01                                                  OBSDEC  =   3.126777777778E+01                                                  OBSGEO-X=   2.225142180269E+06                                                  OBSGEO-Y=  -5.440307370349E+06                                                  OBSGEO-Z=  -2.481029851874E+06                                                  DATE    = '2020-03-03T15:25:51.341000' /Date FITS file was written              ORIGIN  = 'CASA 5.4.0-68'                                                       \n"
