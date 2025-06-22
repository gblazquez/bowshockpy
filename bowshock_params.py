"""
Use this file to define all the the parameters needed to run bowshockpy.

For more information about the physical meaning of some of these parameters, see
the documentation
"""

"""
MODEL OUTPUTS
"""
# Name of the model folder
modelname = f"example_250622"

# Plot 2D bowshock model [True/False]
bs2Dplot = True

# List of the output cube quantites and operations performed onto the cubes.
# The string should follow this format: {quantity}_{operations}
# Available quantities:
#     m: Mass
#     NCO: CO column density
#     tau: Opacity
#     I: Intensity
#     Ithin: Intensity taking into account the optically thin approximation.
#
# Available operations:
#     s: Add source at the origin (just for spatial reference purposes)
#     r: Rotate the cube
#     n: Add gaussian noise
#     c: Convolve with the specified beam
# Operations can be combined, and will be performed following the order of the
# symbols in the string (from left to right).
# Some examples of the strings that could be included in the list outcubes are:
#    "m": Compute the masses in every pixel and channel.
#    "tau_c": Compute the opacities in every pixel and channel, and convolves.
#    "I_rnc": Compute the intensities in every pixel and channel, rotate, add
#    noise, and convolve.
# The list can be left empty if no output cube is desired
# Example of outcubes:
# outcubes = ["m", "m_r", "I_c", "I_nc", "tau_rc", "NCO_rc", "Ithin_rc"]
outcubes = ["I_nc"]

# List of the cubes to which the position-velocity diagrams and moments (0, 1,
# 2, and peak intensity) are going to be performed.
# Example of momentsandpv:
# momentsandpv = ["I_rc", "Ithin_rc"]
momentsandpv = ["I_nc"]

# Verbose messages about the computation? [True/False]
verbose = True

"""
OBSERVER PARAMETERS
"""

# Source distance to the observer [pc]
distpc = 300

# Systemic velocity of the source [km/s]
vsys = + 5

# Source coordinates [deg, deg]
ra_source_deg, dec_source_deg = 51.41198333, 30.73479833


"""
BOWSHOCK PARAMETERS
"""

# Number of bowshocks to model
nbowshocks = 1

# Excitation temperature [K]
Tex = 100

# Background temperature [K]
Tbg = 2.7

# Mean molecular mass per H molecule
muH2 = 2.8

# CO abundance
XCO = 8.5 * 10**(-5)

# The individual bowshock parameters must end in _{bowshock_number}. For example, the jet
# velocity for the third bowshock is vj_3

"""
bowshock 1 [blue]
"""

# Jet inclination angle with respect to the line of sight. If i>90, the jet is
# redshifted, if i<90, it will be blueshifted. [degrees]
i_1 = 180-45

# Jet radius. Set this parameter to zero, the channel maps generator
# are not yet generalized for jet radius>0 [arcsec]
rj_1 = 0

# Characteristic length scale [arcsec]
L0_1 = 0.7

# Distance between the working surface and the source [arcsec]
zj_1 = 3.5

# Jet velocity
vj_1 = 73

# Ambient (or wind) velocity [km/s]
vw_1 = 0

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
SPECTRAL CUBE PARAMETERS
"""

# Number of points to model
nzs = 500

# Number of azimuthal angle phi to calculate the bowshock solution
nphis = 500

# Number of spectral channel maps
nc = 50

# Central velocity of the first channel map [km/s]
vch0 = 40

# Central velocity of the last channel map [km/s]
vchf = +100

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
refpix = [100, 0]

# Angle to rotate the image [degrees]
parot = 0

# Add source to the image at the reference pixel? [True/False]
add_source = False

# Map noise
# Standard deviation of the noise of the map relative to the maximum pixel in the cube, before convolving the cube. The actual noise will be computed after convolving.
maxcube2noise = 15



"""
MOMENTS AND PV PARAMETERS
"""

# Do you want to save the moments and the pv in fits file? [True/False]
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
