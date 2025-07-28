
html-theme.sidebar_secondary.remove: true

=======================
Examples of input file
=======================

Here, some examples of input files are presented. See :doc:`usage<../user_guide/usage>` section for an explanation of how to use this input file, and :doc:`input parameters<../user_guide/inputparams>` secton for a description of the parameters.

You can copy and paste these examples to your local machine, download them from the `examples <https://github.com/gblazquez/bowshockpy/tree/main/examples>`_ folder available in the GitHub repository, or print them directy to your working directory

.. code-block:: console

  (.venv) $ bowshockpy --print-example 1

This will print example 1 to your working directory. Then, you can modify the example file according to your needs. 


Example 1: A redshifted bowshock
---------------------------------------------

This example of input file generates one redshifted bowshock. As specified by the parameter outcube, the output will be a cube of the intensities, with gaussian noise, and convolved (its filename is 'I_nc.fits'). The moment images and the PV diagram along the jet axis will also be computed. The opacity, masses, and CO column densities will also be saved.

.. code-block:: python
  
  https://bowshockpy.readthedocs.io/en/latest/user_guide/inputparams.html
  
  """
  
  """
  MODEL OUTPUTS
  """
  modelname = f"example1"
  
  bs2Dplot = True
  
  outcubes = {
      "intensity": ["add_noise", "convolve", "moments_and_pv"],
      "opacity": [],
      "CO_column_density": [],
      "mass": [],
      }
  
  verbose = True
  
  """
  OBSERVER PARAMETERS
  """
  
  distpc = 400
  
  vsys = + 5
  
  ra_source_deg, dec_source_deg = 84.095, -6.7675
  
  
  """
  BOWSHOCKS PARAMETERS
  """
  
  nbowshocks = 1
  
  Tex = 100
  
  Tbg = 2.7
  
  muH2 = 2.8
  
  J = 3
  
  XCO = 8.5 * 10**(-5)
  
  
  
  """
  bowshock 1 [redshifted]
  """
  
  i_1 = 135
  
  L0_1 = 0.7
  
  zj_1 = 3.5
  
  vj_1 = 73
  
  va_1 = 0
  
  v0_1 = 5
  
  rbf_obs_1 = 1
  
  mass_1 = 0.00015
  
  pa_1 = -20
  
  
  """
  SPECTRAL CUBE PARAMETERS
  """
  
  nzs = 1000
  
  nphis = 500
  
  nc = 50
  
  vch0 = 35
  
  vchf = 65
  
  nxs = 200
  
  nys = 200
  
  xpmax = 4
  
  papv = pa_1
  
  bmaj, bmin = (0.420, 0.287)
  
  pabeam = -17.2
  
  vt = "2xchannel"
  
  tolfactor_vt = 3
  
  CIC = True
  
  refpix = [80, 30]
  
  coordcube = "sky"
  
  parot = 0
  
  sigma_beforeconv = 0.05
  
  maxcube2noise = 0.07
  
  
  """
  MOMENTS AND PV PARAMETERS
  """
  
  savefits = True
  
  saveplot = True
  
  mom1clipping = "5xsigma"
  
  mom2clipping = "4xsigma"
  
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

Example 2: A blueshifted bowshock
---------------------------------------------

This example of input file generates one blueshifted bowshock. As defined by outcube parameter, the intensities will be computed with and without taking into account the optically thin approximation, gaussian noise will be added and the cubes will be convolved. Moments images and the PV diagram along the jet axis will be computed.

.. code-block:: python
  
  https://bowshockpy.readthedocs.io/en/latest/user_guide/inputparams.html
  
  """
  
  """
  MODEL OUTPUTS
  """
  modelname = f"example2"
  
  bs2Dplot = True
  
  outcubes = {
      "intensity": ["add_noise", "convolve", "moments_and_pv"],
      "intensity_opthin": ["add_noise", "convolve", "moments_and_pv"],
      "opacity": [],
      "mass": [],
      }
  
  verbose = True
  
  """
  OBSERVER PARAMETERS
  """
  
  distpc = 400
  
  vsys = + 5
  
  ra_source_deg, dec_source_deg = 84.095, -6.7675
  
  
  """
  BOWSHOCKS PARAMETERS
  """
  
  nbowshocks = 1
  
  Tex = 100
  
  Tbg = 2.7
  
  muH2 = 2.8
  
  J = 3
  
  XCO = 8.5 * 10**(-5)
  
  
  
  """
  bowshock 1 [redshifted]
  """
  
  i_1 = 25
  
  L0_1 = 0.8
  
  zj_1 = 3.5
  
  vj_1 = 80
  
  va_1 = 0
  
  v0_1 = 10
  
  rbf_obs_1 = 1.1
  
  mass_1 = 0.00015
  
  pa_1 = +40
  
  
  """
  SPECTRAL CUBE PARAMETERS
  """
  
  nzs = 1000
  
  nphis = 500
  
  nc = 50
  
  vch0 = -25
  
  vchf = -80
  
  nxs = 200
  
  nys = 200
  
  xpmax = 4
  
  papv = pa_1
  
  bmaj, bmin = (0.420, 0.287)
  
  pabeam = -17.2
  
  vt = "2xchannel"
  
  tolfactor_vt = 3
  
  CIC = True
  
  refpix = [125, 75]
  
  coordcube = "sky"
  
  parot = 0
  
  
  sigma_beforeconv = 0.03
  
  maxcube2noise = 0.07
  
  
  """
  MOMENTS AND PV PARAMETERS
  """
  
  savefits = True
  
  saveplot = True
  
  mom1clipping = "5xsigma"
  
  mom2clipping = "4xsigma"
  
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

Example 3: A side-on bowshock
---------------------------------------------

This example of input file generates a bowhsock that is side-on; that is, in nearly contain in the plane-of-sky and, consequently, has blue- and red-shifted parts. As specified in outcube parameter, the intensities will be convolved and gaussian noise will be added. Also, the moments and the position velocity diagram will be computed. The cubes of the opcities, CO_column densities and masses are going also to be saved.

.. code-block:: python
  
  https://bowshockpy.readthedocs.io/en/latest/user_guide/inputparams.html
  
  """
  
  """
  MODEL OUTPUTS
  """
  modelname = f"example3"
  
  bs2Dplot = True
  
  outcubes = {
      "intensity": ["add_noise", "convolve", "moments_and_pv"],
      "opacity": [],
      "CO_column_density": [],
      "mass": [],
      }
  
  verbose = True
  
  """
  OBSERVER PARAMETERS
  """
  
  distpc = 400
  
  vsys = + 5
  
  ra_source_deg, dec_source_deg = 84.095, -6.7675
  
  
  """
  BOWSHOCKS PARAMETERS
  """
  
  nbowshocks = 1
  
  Tex = 100
  
  Tbg = 2.7
  
  muH2 = 2.8
  
  J = 3
  
  XCO = 8.5 * 10**(-5)
  
  
  
  """
  bowshock 1 [redshifted]
  """
  
  i_1 = 95
  
  L0_1 = 0.7
  
  zj_1 = 3.25
  
  vj_1 = 60
  
  va_1 = 0
  
  v0_1 = 5
  
  rbf_obs_1 = 1
  
  mass_1 = 0.00015
  
  pa_1 = 0
  
  
  """
  SPECTRAL CUBE PARAMETERS
  """
  
  nzs = 1000
  
  nphis = 500
  
  nc = 50
  
  vch0 = 2
  
  vchf = 18
  
  nxs = 200
  
  nys = 200
  
  xpmax = 4.5
  
  papv = pa_1
  
  bmaj, bmin = (0.420, 0.287)
  
  pabeam = -17.2
  
  vt = "2xchannel"
  
  tolfactor_vt = 3
  
  CIC = True
  
  refpix = [100, 0]
  
  coordcube = "sky"
  
  parot = 0
  
  sigma_beforeconv = 0.15
  
  maxcube2noise = None
  
  
  """
  MOMENTS AND PV PARAMETERS
  """
  
  savefits = True
  
  saveplot = True
  
  mom1clipping = "5xsigma"
  
  mom2clipping = "4xsigma"
  
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

Example 4: Several bowshocks in one cube
---------------------------------------------

This example of input file generates two redshifted bowshocks in the same cube. Gaussian noise will be added to the intensity cube and then it will be convolved.  Also, the moments and the position velocity diagram will be computed. The cubes of the opcities and masses are going to be saved also.

.. code-block:: python
  
  https://bowshockpy.readthedocs.io/en/latest/user_guide/inputparams.html
  
  """
  
  """
  MODEL OUTPUTS
  """
  modelname = f"example4"
  
  bs2Dplot = True
  
  outcubes = {
      "intensity": ["add_noise", "convolve", "moments_and_pv"],
      "opacity": [],
      "mass": [],
      }
  
  verbose = True
  
  """
  OBSERVER PARAMETERS
  """
  
  distpc = 400
  
  vsys = + 5
  
  ra_source_deg, dec_source_deg = 84.095, -6.7675
  
  
  """
  BOWSHOCKS PARAMETERS
  """
  
  nbowshocks = 2
  
  Tex = 100
  
  Tbg = 2.7
  
  muH2 = 2.8
  
  J = 3
  
  XCO = 8.5 * 10**(-5)
  
  
  
  """
  bowshock 1 [redshifted]
  """
  
  i_1 = 125
  
  L0_1 = 0.7
  
  zj_1 = 3
  
  vj_1 = 73
  
  va_1 = 0
  
  v0_1 = 4
  
  rbf_obs_1 = 1
  
  mass_1 = 0.00015
  
  pa_1 = -20
  
  """
  bowshock 1 [redshifted]
  """
  
  i_2 = 125
  
  L0_2 = 0.8
  
  zj_2 = 4
  
  vj_2 = 77
  
  va_2 = 0
  
  v0_2 = 4
  
  rbf_obs_2 = 1
  
  mass_2 = 0.0002
  
  pa_2 = -20
  
  
  
  """
  SPECTRAL CUBE PARAMETERS
  """
  
  nzs = 1000
  
  nphis = 500
  
  nc = 50
  
  vch0 = 30
  
  vchf = 57
  
  nxs = 200
  
  nys = 200
  
  xpmax = 5
  
  papv = pa_1
  
  bmaj, bmin = (0.420, 0.287)
  
  pabeam = -17.2
  
  vt = "2xchannel"
  
  tolfactor_vt = 3
  
  CIC = True
  
  refpix = [80, 30]
  
  coordcube = "sky"
  
  parot = 0
  
  sigma_beforeconv = 0.05
  
  maxcube2noise = 0.07
  
  
  """
  MOMENTS AND PV PARAMETERS
  """
  
  savefits = True
  
  saveplot = True
  
  mom1clipping = "5xsigma"
  
  mom2clipping = "4xsigma"
  
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
