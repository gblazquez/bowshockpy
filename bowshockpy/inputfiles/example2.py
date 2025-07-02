"""
Use this input file to define all the the parameters needed to run bowshockpy:

(env.)$ bowshockpy --read <inputfile.py>

For more information about the meaning of some of these parameters, see the documentation: https://bowshockpy.readthedocs.io/en/latest/
"""

"""
MODEL OUTPUTS
"""
# Name of the model folder
modelname = f"example2"

# Plot 2D bowshock model [True/False]
bs2Dplot = True

# Dictionary of the desired output spectral cubes and the operations performed
# over them. The keys of the dictionary are strings indicating the quantities of
# the desired cubes. These are the available quantities of the spectral cubes:
#    "mass": Total mass of molecular hydrogen in solar mass
#    "CO_column_density": Column density of the CO in cm-2.
#    "intensity": Intensity in Jy/beam.
#    "intensity_opthin": Intensity in Jy/beam, using the optically thin approximation.
#    "tau": Opacities
#
# The values of the dictionary are lists of strings indicating the operations to be
# performed over the cube. These are the available operations:
#     "add_source": Add a source at the reference pixel, just for spatial
#     reference purposes.
#     "rotate": Rotate the whole spectral cube by an angle given by parot parameter.
#     "add_noise": Add gaussian noise, defined by maxcube2noise parameter.
#     "convolve": Convolve with a gaussian defined by the parameters bmaj, bmin,
#     and pabeam.
#     "moments_and_pv": Computes the moments 0, 1, and 2, the maximum intensity
#     and the PV diagram.
# The operations will be performed folowing the order of the strings in the list
# (from left to right). The list can be left empty if no operations are desired.
# Examples of outcubes dictionaries:
#
# - The next dictionary will produce 2 cubes, one of the intensities and another with
# the intensities computed with the optically thin approximation. Gaussian
# noise will be added to both of them, then they will be convolved and the
# moments and PV will be computed:
#   outcubes = {
#     "intensity": ["add_noise", "convolve", "moments_and_pv"],
#     "intensity_opthin": ["add_noise", "convolve", "moments_and_pv"],
#     }
#
# - The next dictionary will produce 4 cubes: one with the intensities with
# noise, convolved and with moments and pv computed, another with only the
# intensities (without operations applied to it), and two more cubes with the
# masses and opacities.
#   outcubes = {
#     "intensity": ["add_noise", "convolve", "moments_and_pv"],
#     "intensity": [],
#     "mass": [],
#     "opacity": [],
#     }
#
outcubes = {
    "intensity": ["add_noise", "convolve", "moments_and_pv"],
    "intensity_opthin": ["add_noise", "convolve", "moments_and_pv"],
    "opacity": ["convolve"],
    "opacity": [],
    "CO_column_density": ["convolve"],
    "mass": [],
    }

# Verbose messages about the computation? [True/False]
verbose = True

"""
OBSERVER PARAMETERS
"""

# Source distance to the observer [pc]
distpc = 300

# Systemic velocity of the source [km/s]
vsys = + 0

# Source coordinates [deg, deg]
ra_source_deg, dec_source_deg = 51.41198333, 30.73479833


"""
BOWSHOCK PARAMETERS
"""

# Number of bowshocks to model
nbowshocks = 2

# Excitation temperature [K]
Tex = 100

# Background temperature [K]
Tbg = 2.7

# Mean molecular mass per H molecule
muH2 = 2.8

# CO rovibrational line
J = "3-2"

# CO abundance
XCO = 8.5 * 10**(-5)

# The individual bowshock parameters must end in _{bowshock_number}. For example, the jet
# velocity for the third bowshock is vj_3

"""
bowshock 1 [redshifted]
"""

# Jet inclination angle with respect to the line of sight. If i>90, the jet is
# redshifted, if i<90, it will be blueshifted. [degrees]
i_1 = 180-45

# Characteristic length scale [arcsec]
L0_1 = 0.7

# Distance between the working surface and the source [arcsec]
zj_1 = 3.5

# Jet velocity [km/s]
vj_1 = 73

# Ambient (or surrounding wind) velocity [km/s]
va_1 = 0

# Velocity at which the material is ejected from the internal working surface [km/s]
v0_1 = 5

# Final radius of the bowshock [arcsec]. Set None if you want to end the
# bowshock model at the theoretical final radius (see eq. 11 from Tabone et al.
# 2018)
rbf_obs_1 = 1

# Total mass of the bowshock [solar masses]
mass_1 = 0.00031 * 1.5

# Position angle [deg]
pa_1 = -20

"""
bowshock 2 [redshifted]
"""

# Jet inclination angle with respect to the line of sight. If i>90, the jet is
# redshifted, if i<90, it will be blueshifted. [degrees]
i_2 = 180-45

# Characteristic length scale [arcsec]
L0_2 = 0.8

# Distance between the working surface and the source [arcsec]
zj_2 = 4.5

# Jet velocity [km/s]
vj_2 = 80

# Ambient (or surrounding wind) velocity [km/s]
va_2 = 0

# Velocity at which the material is ejected from the internal working surface [km/s]
v0_2 = 7

# Final radius of the bowshock [arcsec]. Set None if you want to end the
# bowshock model at the theoretical final radius (see eq. 11 from Tabone et al.
# 2018)
rbf_obs_2 = 1

# Total mass of the bowshock [solar masses]
mass_2 = 0.00025 

# Position angle [deg]
pa_2 = -20




"""
SPECTRAL CUBE PARAMETERS
"""

# Number of points to model
nzs = 100

# Number of azimuthal angle phi to calculate the bowshock solution
nphis = 500

# Number of spectral channel maps
nc = 50

# Central velocity of the first channel map [km/s]
vch0 = 30

# Central velocity of the last channel map [km/s]
vchf = 70

# Number of pixels in the x and y axes
nxs, nys = (200, 200)

# Physical size of the channel maps along the x axis [arcsec]
xpmax = 5

# Position angle used to calculate the PV [degrees]
papv = pa_1

# Beam size [arcsec]
bmaj, bmin = (0.420, 0.287)

# Beam position angle [degrees]
pabeam = -17.2

# Thermal+turbulent line-of-sight velocity dispersion [km/s] If
# thermal+turbulent line-of-sight velocity dispersion is smaller than the
# instrumental spectral resolution, vt should be the spectral resolution.
# It can be also set to a integer times the channel width (e.g., "2xchannel")
vt = "2xchannel"

# Cloud in Cell interpolation? [True/False]
CIC = True

# Neighbour channel maps around a given channel map with vch will stop being
# populated when their difference in velocity with respect to vch is higher than
# this factor times vt. The lower the factor, the quicker will be the code, but
# the total mass will be underestimated. If vt is not None, compare the total
# mass of the output cube with the 'mass' parameter that the user has defined
tolfactor_vt = 3

# Reference pixel [[int, int] or None]
# Pixel coordinates (zero-based) of the source, i.e., the origin from which the
# distances are measured. The first index is the R.A. axis, the second is the
# Dec. axis.
refpix = [80, 30]

# Angle to rotate the image [degrees]
parot = 0

# Map noise
# Standard deviation of the noise of the map, before convolution. Set to None if maxcube2noise is used.
sigma_beforeconv = 0.1

# Standard deviation of the noise of the map, before convolution, relative to
# the maximum pixel in the cube. The actual noise will be computed after
# convolving. This parameter would not be used if sigma_beforeconve is not
# None.
maxcube2noise = 0.07



"""
MOMENTS AND PV PARAMETERS
"""

# Do you want to save the moments and the pv in fits format? [True/False]
savefits = True

# Do you want to save a figure of the moments and the PV? [True/False]
saveplot = True

# Clipping for moment 1.
mom1clipping = "5xsigma"

# Clipping for moment 2.
mom2clipping = "4xsigma"

# Set the maximum, central, and minimum value to show in the plot of the moments
# and pv-diagram along the jet axis
mom0values = {
    "vmax": None,
    "vcenter": None,
    "vmin": None,
}

mom1values = {
    "vmax": None,
    "vcenter": None,
    "vmin": None,
}

mom2values = {
    "vmax": None,
    "vcenter": None,
    "vmin": None,
}

mom8values = {
    "vmax": None,
    "vcenter": None,
    "vmin": None,
}

pvvalues = {
    "vmax": None,
    "vcenter": None,
    "vmin": None,
}
