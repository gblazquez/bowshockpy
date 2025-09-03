"""
Use this input file to define all the the parameters needed to run bowshockpy:

$ bowshockpy --read <inputfile.py>

For more information about the meaning of these parameters, see the
documentation:

https://bowshockpy.readthedocs.io/en/latest/user_guide/inputparams.html

"""

"""
MODEL OUTPUTS
"""
# Folder name where the outputs of the modellling are going to be stored. If
# it does not exist, it will be created automatically.
modelname = f"example5"

# Dictionary indicating the desired output spectral cubes and the operations
# performed over them. The keys of the dictionary are strings indicating the
# quantities of the desired cubes. These are the available quantities of the
# spectral cubes:
#
#  - "mass": Total mass of molecular hydrogen in solar mass
#  - "Ntot_column_density": Total (H2 + heavier components) column density in
#  cm-2.
#  - "mol_column_density": Column density of the emitting molecule in cm-2.
#  - "intensity": Intensity in Jy/beam.
#  - "tau": Opacities.
#
# The values of the dictionary are lists of strings indicating the operations to
# be performed over the cube. These are the available operations:
#
#  - "add_source": Add a source at the reference pixel, just for spatial
#  reference purposes.
#  - "rotate": Rotate the whole spectral cube by an angle
#  given by parot parameter.
#  - "add_noise": Add Gaussian noise, defined by maxcube2noise parameter.
#  - "convolve": Convolve with a Gaussian defined by the parameters bmaj, bmin,
#  and pabeam.
#  - "moments_and_pv": Computes the moments 0, 1, and 2, the maximum intensity
#  and the PV diagram.
#
# The operations will be performed folowing the order of the strings in the list
# (from left to right). The list can be left empty if no operations are desired.
#
# For example, the following dictionary for the outcubes parameter,
#
# outcubes = {
#     "intensity": ["add_noise", "convolve", "moments_and_pv"],
#     "opacity": [],
#     "mol_column_density": ["convolve"],
#     "mass": [],
# }
# will save 4 spectral cubes in fits format. The first one are the intensities
# with Gaussian noise added, it will be convolved, and the moments and PV
# diagrams will be computed; the second cube will be the opacity; the third will
# be the mol_column_density, which will be convolved; and the forth cube will be
# the masses.
outcubes = {
    "intensity": ["convolve", "moments_and_pv"],
    "opacity": [],
    "mol_column_density": [],
    "mass": [],
    }

# Set True to verbose messages about the computation.
verbose = True

"""
OBSERVER PARAMETERS
"""

# Source distance to the observer [pc]
distpc = 400

# Systemic velocity of the source [km/s]
vsys = + 5

# Source coordinates [deg, deg]
ra_source_deg, dec_source_deg = 84.095, -6.7675


"""
BOWSHOCKS PARAMETERS
"""

# Number of bowshocks to model
nbowshocks = 1

# Upper level of the rotational transition (e.g. 3 for the "J=3->2" transition).
J = 4

# Frequency of the transition [GHz]
nu = 576.2679305

# Abundance relative to the molecular hydrogen.
abund = 4 * 10**(-5)

# Mean molecular mass per hydrogen molecule
meanmolmass = 2.8

# Permanent dipole moment [Debye]
mu = 0.112

# Excitation temperature [K]
Tex = 100

# Background temperature [K]
Tbg = 2.7


# You can model several bowshocks in the same spectral cube. The number of
# bowshocks are given by **nbowshocks** parameter. The following parameters
# should be defined for each bowshock, subtituting "n" with the bowshock index
# (e.g., if 4 bowshocks are included in the model, one should define **vj_1**,
# **vj_2**, **vj_3**, and **vj_4**, and similarly with the rest of parameters).

"""
bowshock 1 [redshifted]
"""

# Inclination angle of the bowshock symmetry axis with respect to the line of
# sight. If i>90, the bowshock is redshifted, if i<90, it will be blueshifted
# [degrees].
i_1 = 135

# Characteristic length scale [arcsec]
L0_1 = 0.7

# Distance between the internal working surface and the source [arcsec]
zj_1 = 3.5

# Jet velocity [km/s]
vj_1 = 73

# Ambient (or surrounding wind) velocity [km/s]
va_1 = 0

# Velocity at which the material is ejected sideways from the internal working
# surface [km/s]
v0_1 = 5

# Final radius of the bowshock [arcsec]. Set None if you want to end the
# bowshock model at the theoretical final radius (see eq. 11 from Tabone et
# al. 2018).
rbf_obs_1 = 1

# Total mass of the bowshock [solar masses]
mass_1 = 0.00015

# Position angle [deg]
pa_1 = -20


"""
SPECTRAL CUBE PARAMETERS
"""

# Number of points to model along the direction of the symmetry axis (z-axis).
nzs = 1000

# Number of azimuthal angles to calculate the bowshock solution at each model
# point in the z-axis.
nphis = 500

# Number of spectral channel maps.
nc = 50

# Central velocity of the first channel map [km/s].
vch0 = 35

# Central velocity of the last channel map [km/s]. Set to None if chanwidth is
# used.
vchf = 65

# Width of the velocity channel [km/s]. If chanwidth>0, then vch0<vchf, if
# chanwidth<0, then vch0>vchf. Set to None if vchf is used.
chanwidth = None

# Number of pixels in the right ascension axis.
nxs = 200

# Number of pixels in the declination axis.
nys = 200

# Physical size of the channel maps along the right ascension axis [arcsec].
xpmax = 4

# Position angle used to calculate the PV [degrees]
papv = pa_1

# Beam size [arcsec]
bmaj, bmin = (0.420, 0.287)

# Beam position angle [degrees]
pabeam = -17.2

# Thermal+turbulent line-of-sight velocity dispersion [km/s] If
# thermal+turbulent line-of-sight velocity dispersion is smaller than the
# instrumental spectral resolution, vt should be the spectral resolution.
# It can be also set to a integer times the channel width (e.g., "2xchannel").
vt = "2xchannel"

# The masses corresponding to a channel map are spread along the whole cube in
# the velocity axis following a Gaussian distribution, being vt parameter the
# standard deviation of the Gaussian. tolfactor_vt parameter truncates the
# Gaussian distribution at vt * tolfactor_vt in order to make the computation
# substatially faster. A low tolfactor_vt can result in a warning reporting an
# underestimation of the total mass of the model.
tolfactor_vt = 3

# Set to True to perform 2D Cloud in Cell interpolation along the spatial
# dimensions. If False, a Nearest Grid Point method will be perform.
cic = True

# Pixel coordinates (zero-based) of the source, i.e., the origin from which
# the distances are measured. The first index is the right ascension axis,
# the second is the declination axis [[int, int] or None].
refpix = [80, 30]

# Set to "sky" in order to set the cube headers in sky coordinates, or "offset"
# if you prefer them in offsets relative to the origin (the source).
coordcube = "sky"

# Angle to rotate the image [degrees]
parot = 0

# Standard deviation of the noise of the map, before convolution. Set to None if
# maxcube2noise is used.
sigma_beforeconv = 0.05

# Standard deviation of the noise of the map, before convolution, relative to
# the maximum pixel in the cube. The actual noise will be computed after
# convolving. This parameter would not be used if sigma_beforeconve is not None.
maxcube2noise = 0.07


"""
MOMENTS AND PV PARAMETERS
"""

# Set to True in order save the moments and the PV in fits format.
savefits = True

# Set to True in order to save a figure of the moments and the PV [True/False].
saveplot = True

# Clipping for moment 1 as a function of the standard deviation of noise in the
# image (e.g., "5xsigma").
mom1clipping = "5xsigma"

# Clipping for moment 2 as a function of the standard deviation of noise in the
# image (e.g., "4xsigma").
mom2clipping = "4xsigma"

# Dictionary with the maximum, central, and minimum value to show in the plot of
# the moment 0. If the dictionary value is None for vmax, vcenter, or vmin, then
# the maximum, central, or the minimum value of the moment image will be
# considered, respectively. Example: mom0values = {"vmax": None, "vcenter":
# None, "vmin": 0,}.
mom0values = {
    "vmax": None,
    "vcenter": None,
    "vmin": None,
}

# Dictionary with the maximum, central, and minimum value to show in the plot of
# the moment 1. If the dictionary value is None for vmax, vcenter, or vmin, then
# the maximum, central, or the minimum value of the moment image will be
# considered, respectively. Example: mom1values = {"vmax": 60, "vcenter": 20,
# "vmin": 0,}.
mom1values = {
    "vmax": None,
    "vcenter": None,
    "vmin": None,
}

# Dictionary with the maximum, central, and minimum value to show in the plot of
# the moment 2. If the dictionary value is None for vmax, vcenter, or vmin, then
# the maximum, central, or the minimum value of the moment image will be
# considered, respectively. Example: mom2values = {"vmax": None, "vcenter":
# None, "vmin": None,}.
mom2values = {
    "vmax": None,
    "vcenter": None,
    "vmin": None,
}

# Dictionary with the maximum, central, and minimum value to show in the plot
# of the maximum value along the velocity axis. If the dictionary value is
# None for vmax, vcenter, or vmin, then the maximum, central, or the minimum
# value of the moment image will be considered, respectively. Example:
# maxintensvalues = {"vmax": None, "vcenter": None, "vmin": None,}.
maxintensvalues = {
    "vmax": None,
    "vcenter": None,
    "vmin": None,
}

# Set the maximum, central, and minimum value to show in the plot of the moments
# and PV-diagram along the jet axis. If the dictionary value is None for vmax,
# vcenter, or vmin, then the maximum, central, or the minimum value of the
# position velocity diagram will be considered, respectively. Example: pvvalues
# = {"vmax": None, "vcenter": None, "vmin": None,}.
pvvalues = {
    "vmax": None,
    "vcenter": None,
    "vmin": None,
}


"""
CUSTOM TRANSITION MODEL AND RADIATIVE TRANSFER
(Optional)
"""
# BowshockPy, by default, is able to compute the intensities of a low-J
# rotational transition of a linear molecule, neglecting the population of
# excitated vibrational levels and centrifugal distortions of the molecule.
# Also, Local Thermodynamic Equilibrium. If you need a custom model the
# molecular transition or a different radiative transfer, you can implement
# them in the in the following lines.
#
# First, import the necessary modules of your script
# import astropy.units as u
# import astropy.constants as const
#
# You can define as many functions or classes as you need to implement your
# model. Nonetheless, BowshockPy will need only two:
#
# Custom function to compute the opacities. Should be a function only of the
# column densities per velocity bin
# def tau_custom_function(dNmoldv):
#     """
#     Custom function to compute the opacities from the column densities per
#     velocity bin
#
#     Parameters
#     ----------
#     dNmoldv : astropy.units.Quantity
#         Column density per velocity bin
#
#     Returns
#     -------
#     tau : float
#         Opacity
#     """
#     pass
#
# Custom function to compute the intensities from the opacities
# def Inu_custom_function(tau):
#     """
#     Custom function to perform the radiative transfer from the opacities
#
#     Parameters
#     ----------
#     tau : float
#         Opacity
#
#     Returns
#     -------
#     Inu : astropy.units.Quantity
#         Intensity
#     """
#     pass


# Import the necessary function of your script
import bowshockpy.radtrans as rt
import astropy.constants as const
import astropy.units as u


# Define the functions that you think are useful for your modelling
def Ej(j, B0, D0):
    """
    Energy state of a rotational transition of a linear molecule, taking
    into account the first order centrifugal distortion

    Parameters
    ----------
    j : int
        Rotational level
    B0 : astropy.units.quantity
        Rotation constant
    D0 : astropy.units.quantity
        First order centrifugal distortion constant

    Returns
    -------
    astropy.units.quantity
        Energy state of a rotator
    """
    return const.h * (B0 * j * (j+1) - D0 * j**2 * (j+1)**2)

def gj(j):
    """
    Degeneracy of the level j at which the measurement was made. For a
    linear molecule, g = 2j + 1

    Parameters
    ----------
    j : int
        Rotational level

    Returns
    -------
    int
        Degeneracy of the level j
    """
    return 2*j + 1

def muj_jm1(j, mu_dipole):
    """
    Computes the dipole moment matrix element squared for rotational
    transition j->j-1

    Parameters
    ----------
    j : int
        Rotational level
    mu_dipole : astropy.units.quantity
        Permanent dipole moment of the molecule
    """
    return mu_dipole * (j / (2*j + 1))**0.5

# Only the next two functions will be read by BowshockPy

def tau_custom_function(dNmoldv):
    """
    Custom function to compute the opacities from the column densities per
    velocity bin

    Parameters
    ----------
    dNmoldv : astropy.units.Quantity
        Column density per velocity bin

    Returns
    -------
    tau : float
        Opacity
    """

    B0 = 57.62 * u.GHz # nu / (2J)
    D0 = B0 * 2 * 10**(-5)
    mu_ul = muj_jm1(J, mu*u.Debye)
    # We can perform the calculation of the partition function and tau from the
    # scratch, or we can use the function tau_func from bowshockpy.radtrans
    # module, which computes internally the partition function from the
    # user defined function Ei(i, *args), which computes the energy of level i.
    tau = rt.tau_func(
        dNmoldv=dNmoldv,
        nu=nu*u.GHz,
        Tex=Tex*u.K,
        i=J,
        Ei=Ej,
        gi=gj,
        mu_ul=mu_ul,
        Ei_args=(B0, D0), # pass all the extra arguments to Ei
        gi_args=(),
    )
    return tau


def Inu_custom_function(tau):
    """
    Computes the intensity through the radiative transfer equation. We assume
    optically thin emission

    Parameters
    ----------
    tau : float
        Opacity

    Returns
    -------
    astropy.units.quantity
        Intensity (energy per unit of area, time, frequency and solid angle)
    """
    Inu = rt.Bnu_func(nu*u.GHz, Tex*u.K) * tau
    return Inu
