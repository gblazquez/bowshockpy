

=======================
Examples of input file
=======================

The quickest way to use ``BowshockPy`` is to run it from the terminal and specify an input file containing all the parameters needed to generate the models. You can find here four different examples of input files:

- `Example 1: A redshifted bowshock <Example 1: A redshifted bowshock>`_
- `Example 2: A blueshifted bowshock <Example 2: A blueshifted bowshock>`_
- `Example 3: A side-on bowshock <Example 3: A side-on bowshock>`_
- `Example 4: Several bowshocks in one cube <Example 4: Several bowshocks in one cube>`_
- `Example 5: Custom model for the molecular transition <Example 5: Custom model for the molecular transition>`_

You can either copy and paste these examples to your local machine, download them from the `examples <https://github.com/gblazquez/bowshockpy/tree/main/examples>`_ folder available in the `GitHub repository <https://github.com/gblazquez/bowshockpy>`_, or print them directy to your working directory

.. code-block:: console

  $ bowshockpy --print example1.py 

This will copy example 1 (example1.py) to your working directory. Then, you can modify the example file according to your needs and run ``bowshockpy``


.. code-block:: console

  $ bowshockpy --read example1.py


See :doc:`usage<../user_guide/usage>` section for more detailed explanation of how to run bowshockpy, and :doc:`input parameters<../user_guide/inputparams>` section for a description of the parameters.


Example 1: A redshifted bowshock
----------------------------------------------------------------

This example of input file generates one redshifted bowshock. As specified by the parameter outcube, the output will be a cube of the intensities, with Gaussian noise, and convolved (its filename is 'I_nc.fits'). The moment images and the PV diagram along the jet axis will also be computed. The opacity, masses, and CO column densities will also be saved.

.. code-block:: python
  
  """
  MODEL OUTPUTS
  """
  modelname = f"example1"
  
  outcubes = {
      "intensity": ["add_noise", "convolve", "moments_and_pv"],
      "opacity": [],
      "total_column_density": [],
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
  
  J = 3
  
  nu = 345.79598990
  
  abund = 8.5 * 10**(-5)
  
  meanmolmass = 2.8
  
  mu = 0.112
  
  Tex = 100
  
  Tbg = 2.7
  
  
  
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
  
  chanwidth = None
  
  nxs = 200
  
  nys = 200
  
  xpmax = 4
  
  papv = pa_1
  
  bmaj, bmin = (0.420, 0.287)
  
  pabeam = -17.2
  
  vt = "2xchannel"
  
  tolfactor_vt = 3
  
  cic = True
  
  refpix = [80, 30]
  
  coordcube = "sky"
  
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
  
  maxintensvalues = {
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
----------------------------------------------------------------

This example of input file generates one blueshifted bowshock. As defined by outcube parameter, the intensities will be computed with and without taking into account the optically thin approximation, Gaussian noise will be added and the cubes will be convolved. Moments images and the PV diagram along the jet axis will be computed.

.. code-block:: python
  
  """
  MODEL OUTPUTS
  """
  modelname = f"example2"
  
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
  
  nbowshocks = 1
  
  J = 3
  
  nu = 345.79598990
  
  abund = 8.5 * 10**(-5)
  
  meanmolmass = 2.8
  
  mu = 0.112
  
  Tex = 100
  
  Tbg = 2.7
  
  
  
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
  
  chanwidth = None
  
  nxs = 200
  
  nys = 200
  
  xpmax = 4
  
  papv = pa_1
  
  bmaj, bmin = (0.420, 0.287)
  
  pabeam = -17.2
  
  vt = "2xchannel"
  
  tolfactor_vt = 3
  
  cic = True
  
  refpix = [125, 75]
  
  coordcube = "sky"
  
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
  
  maxintensvalues = {
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
----------------------------------------------------------------

This example of input file generates a bowhsock that is side-on; that is, in nearly contain in the plane-of-sky and, consequently, has blue- and red-shifted parts. As specified in outcube parameter, the intensities will be convolved and Gaussian noise will be added. Also, the moments and the position velocity diagram will be computed. The cubes of the opcities, CO_column densities and masses are going also to be saved.

.. code-block:: python
  
  """
  MODEL OUTPUTS
  """
  modelname = f"example3"
  
  outcubes = {
      "intensity": ["add_noise", "convolve", "moments_and_pv"],
      "opacity": [],
      "emitting_molecule_column_density": [],
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
  
  J = 3
  
  nu = 345.79598990
  
  abund = 8.5 * 10**(-5)
  
  meanmolmass = 2.8
  
  mu = 0.112
  
  Tex = 100
  
  Tbg = 2.7
  
  
  
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
  
  chanwidth = None
  
  nxs = 200
  
  nys = 200
  
  xpmax = 4.5
  
  papv = pa_1
  
  bmaj, bmin = (0.420, 0.287)
  
  pabeam = -17.2
  
  vt = "2xchannel"
  
  tolfactor_vt = 3
  
  cic = True
  
  refpix = [100, 0]
  
  coordcube = "sky"
  
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
  
  maxintensvalues = {
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
----------------------------------------------------------------

This example of input file generates two redshifted bowshocks in the same cube. Gaussian noise will be added to the intensity cube and then it will be convolved.  Also, the moments and the position velocity diagram will be computed. The cubes of the opcities and masses are going to be saved also.

.. code-block:: python
  
  """
  MODEL OUTPUTS
  """
  modelname = f"example4"
  
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
  
  J = 3
  
  nu = 345.79598990
  
  abund = 8.5 * 10**(-5)
  
  meanmolmass = 2.8
  
  mu = 0.112
  
  Tex = 100
  
  Tbg = 2.7
  
  
  
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
  
  chanwidth = None
  
  nxs = 200
  
  nys = 200
  
  xpmax = 5
  
  papv = pa_1
  
  bmaj, bmin = (0.420, 0.287)
  
  pabeam = -17.2
  
  vt = "2xchannel"
  
  tolfactor_vt = 3
  
  cic = True
  
  refpix = [80, 30]
  
  coordcube = "sky"
  
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
  
  maxintensvalues = {
      "vmax": None,
      "vcenter": None,
      "vmin": None,
  }
  
  pvvalues = {
      "vmax": None,
      "vcenter": None,
      "vmin": None,
  }

Example 5: Custom model for the molecular transition
----------------------------------------------------------------

This example of input file generates one redshifted. A custom model for the transition is added at the end. As specified by the parameter outcube, the output will be a cube of the intensities, with Gaussian noise, and convolved (its filename is 'I_nc.fits'). The moment images and the PV diagram along the jet axis will also be computed. The opacity, masses, and column densities of the emitting molecule will also be saved.

.. code-block:: python
  
  """
  MODEL OUTPUTS
  """
  modelname = f"example5"
  
  outcubes = {
      "intensity": ["convolve", "moments_and_pv"],
      "opacity": [],
      "emitting_molecule_column_density": [],
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
  
  J = 4
  
  nu = 576.2679305
  
  abund = 4 * 10**(-5)
  
  meanmolmass = 2.8
  
  mu = 0.112
  
  Tex = 100
  
  Tbg = 2.7
  
  
  
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
  
  chanwidth = None
  
  nxs = 200
  
  nys = 200
  
  xpmax = 4
  
  papv = pa_1
  
  bmaj, bmin = (0.420, 0.287)
  
  pabeam = -17.2
  
  vt = "2xchannel"
  
  tolfactor_vt = 3
  
  cic = True
  
  refpix = [80, 30]
  
  coordcube = "sky"
  
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
  
  maxintensvalues = {
      "vmax": None,
      "vcenter": None,
      "vmin": None,
  }
  
  pvvalues = {
      "vmax": None,
      "vcenter": None,
      "vmin": None,
  }
  
  
  """
  CUSTOM TRANSITION MODEL AND RADIATIVE TRANSFER
  (Optional)
  """
  
  
  import bowshockpy.radtrans as rt
  import astropy.constants as const
  import astropy.units as u
  
  
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
