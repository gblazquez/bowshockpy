Input parameters
================

In this section, the input parameters that ``bowshockpy`` needs are described. You can either define these parameters in an input file (the easiest way, see :doc:`input file examples<../examples/examples_inputfile>`), or import ``bowshockpy`` as a python package and define the parameters in a dictionary that would be needed as an input in order instatiate the clases (the most flexibe way, see :doc:`modular usage examples<../examples/example_module>`).


Model output parameters
-----------------------

These parameters define the desired outputs:

*modelname* (str)
    Name of the folder where the outputs of the models are going to be stored. 

*bs2Dplot* (boolean)
    Set to True to plot a 2D bowshock model.

*outcubes* (dict)
    Dictionary indicating the desired output spectral cubes and the operations performed over them. The keys of the dictionary are strings indicating the quantities of the desired cubes. These are the available quantities of the spectral cubes:

    - "mass": Total mass of molecular hydrogen in solar mass
    - "CO_column_density": Column density of the CO in cm-2.
    - "intensity": Intensity in Jy/beam.
    - "intensity_opthin": Intensity in Jy/beam, using the optically thin approximation.
    - "tau": Opacities.

    The values of the dictionary are lists of strings indicating the operations to be performed over the cube. These are the available operations:

    - "add_source": Add a source at the reference pixel, just for spatial reference purposes.
    - "rotate": Rotate the whole spectral cube by an angle given by parot parameter.
    - "add_noise": Add gaussian noise, defined by maxcube2noise parameter.
    - "convolve": Convolve with a gaussian defined by the parameters bmaj, bmin, and pabeam.
    - "moments_and_pv": Computes the moments 0, 1, and 2, the maximum intensity and the PV diagram.

    The operations will be performed folowing the order of the strings in the list (from left to right). The list can be left empty if no operations are desired.
    
    For example, the following dictionary for the outcubes parameter,

    .. code-block:: python
     
        outcubes = {
            "intensity": ["add_noise", "convolve", "moments_and_pv"],
            "opacity": [],
            "CO_column_density": ["convolve"],
            "mass": [],
        }

    will save 4 spectral cubes in fits format. The first one are the intensities with gaussian noise added, it will be convolved, and the moments and PV diagrams will be computed; the second cube will be the opacity; the third will be the CO_column_density, which will be convolved; and the forth cube will be the masses. The first spectral cube will be named I_nc.fits, the second tau.fits, the third NCO_c.fits, and the fourth m.fits. See :doc:`outputs<outputs>` section for a full description of the outputs and the abbreviations used in the filenames of each fits file.


*verbose* (bolean)
    Set True to verbose messages about the computation


Observer parameters
-------------------

These parameters define the observer properties:

*distpc* (float)
    Source distance to the observer [pc]
*vsys* (float)
    Systemic velocity of the source [km/s]
*ra_source_deg* (float)
    Source right ascension [deg]
*dec_source_deg* (float)
    Source declination [deg]


Bowshock parameters
-------------------

These parameters common to all the bowshocks that are going to be generated:

*nbowshocks* (int)
    Number of bowshocks to model
*Tex* (float)
    Excitation temperature [K]
*Tbg* (float)
    Background temperature [K]
*muH2* (float)
    Mean molecular mass per H molecule
*J* (str)
    CO rotational transition (e.g. "3-2")
*XCO* (float)
    CO abundance

``bowhsockpy`` allows to model several bowshocks in the same spectral cube. The number of bowshocks are given by *nbowshocks* parameter. The following paramers should be defined for each bowshock, subtituting "n" with the bowshock index (e.g., if 4 bowshocks are included in the model, one should define vj_1, vj_2, vj_3, and vj_4, and similarly with the rest of parameters).

*i_n* (foat)
    Inclination angle with respect to the line of sight. If i>90, the bowshock is redshifted, if i<90, it will be blueshifted [degrees].
 
*L0_n* (float)
    Characteristic length scale [arcsec]

*zj_n* (float)
    Distance between the working surface and the source [arcsec]

*vj_n* (float)
    Jet velocity [km/s]

*va_n* (float)
    Ambient (or surrounding wind) velocity [km/s]

*v0_n* (float) 
    Velocity at which the material is ejected from the internal working surface [km/s]

*rbf_obs_1* (float)
    Final radius of the bowshock [arcsec]. Set None if you want to end the bowshock model at the theoretical final radius (see eq. 11 from Tabone et al. 2018)

*mass_n* (float)
    Total mass of the bowshock [solar masses]

*pa_n* (float)
    Position angle [deg]


Spectral cube parameters
-------------------------

*nzs* (int)
    Number of points to model

*nphis* (int)
    Number of azimuthal angle phi to calculate the bowshock solution

*nc* (int)
    Number of spectral channel maps

*vch0* (float)
    Central velocity of the first channel map [km/s]

*vchf* (float)
    Central velocity of the last channel map [km/s]

*nxs* (int)
    Number of pixels in the right ascension axis.

*nys* (int)
    Number of pixels in the declination axis. 

*xpmax* (float)
    Physical size of the channel maps along the x axis [arcsec]

*papv* (float)
    Position angle used to calculate the PV [degrees]

*bmaj* (tupple)
    Beam major axis [arcsec]

*bmin* (tupple)
    Beam minor axis [arcsec]

*pabeam* (float)
    Beam position angle [degrees]

*vt* (str or float)
    Thermal+turbulent line-of-sight velocity dispersion [km/s] If thermal+turbulent line-of-sight velocity dispersion is smaller than the instrumental spectral resolution, vt should be the spectral resolution. It can be also set to a integer times the channel width (e.g., "2xchannel")

*tolfactor_vt* (float)
    If the velocity dispersion vt is not 0, the masses corresponding to a channel
    maps are spread along the whole cube in the velocity axis following a Gaussian
    distribution. tolfactor_vt parameter truncates the gaussian distribution at
    vt*tolfactor_vt in order to make the computation substatially faster. A low
    tolfactor_vt can result in a warning reporting an underestimation of the total
    mass of the model.

*CIC* (bolean)
    Set to True to perform 2D Cloud in Cell interpolation along the spatial
    dimensions. If False, a nearest neighbour association will be perform.
    
*refpix* (list or None)
    Pixel coordinates (zero-based) of the source, i.e., the origin from which the distances are measured. The first index is the R.A. axis, the second is the  Dec. axis [[int, int] or None] 

*coordcube* ("sky" or "offset")
    Set to "sky" if you would like to set the cube headers in sky coordinates, or "offset" if you prefer them in offsets relative to the origin (the source).

*parot* (float)
    Angle to rotate the image [degrees]

*sigma_beforeconv* (float)
    Standard deviation of the noise of the map, before convolution. Set to None if maxcube2noise is used.

*maxcube2noise* (float)
    Standard deviation of the noise of the map, before convolution, relative to the maximum pixel in the cube. The actual noise will be computed after convolving. This parameter would not be used if sigma_beforeconve is not None.


Moments and PV parameters
-------------------------

*savefits* (bolean)
    Set to True in order save the moments and the PV in fits file format.

*saveplot* (bolean)
    Do you want to save a figure of the moments and the PV? [True/False]

*mom1clipping* (str)
    Clipping for moment 1 as a function of the standard deviation of noise in the image (e.g., "5xsigma").

*mom2clipping* (str)
    Clipping for moment 2 as a function of the standard deviation of noise in the image (e.g., "4xsigma").

*mom0values* (dict)
    Dictionary with the maximum, central, and minimum value to show in the plot
    of the moment 0. If the dictionary value is None for vmax, vcenter, or vmin,
    then the maximum, central, or the minimum value of the moment image will be
    considered, respectively. Example: mom0values = {"vmax": None, "vcenter": None,
    "vmin": 0,}. 

*mom1values* (dict)
    Dictionary with the maximum, central, and minimum value to show in the plot
    of the moment 1. If the dictionary value is None for vmax, vcenter, or vmin,
    then the maximum, central, or the minimum value of the moment image will be
    considered, respectively. Example: mom1values = {"vmax": 60, "vcenter": 20,
    "vmin": 0,}. 
    
*mom2values* (dict)
    Dictionary with the maximum, central, and minimum value to show in the plot
    of the moment 2. If the dictionary value is None for vmax, vcenter, or vmin,
    then the maximum, central, or the minimum value of the moment image will be
    considered, respectively. Example: mom2values = {"vmax": None, "vcenter": None,
    "vmin": None,}. 

*mom8values* (dict)
    Dictionary with the maximum, central, and minimum value to show in the plot
    of the maximum value along the velocity axis. If the dictionary value is
    None for vmax, vcenter, or vmin, then the maximum, central, or the minimum
    value of the moment image will be considered, respectively. Example:
    mom2values = {"vmax": None, "vcenter": None, "vmin": None,}. 

*pvvalues* (dict) 
    Set the maximum, central, and minimum value to show in the plot of the
    moments and PV-diagram along the jet axis. If the dictionary value is None
    for vmax, vcenter, or vmin, then the maximum, central, or the minimum value
    of the position velocity diagram will be considered, respectively. Example:
    pvvalues = {"vmax": None, "vcenter": None, "vmin": None,}
    
