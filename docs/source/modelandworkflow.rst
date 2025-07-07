How bowhsockpy works
====================

The bowshock model
------------------

``bowshockpy`` computes synthetic spectral cubes, PV diagrams, and moment images for a simple analytical jet-driven bowshock model, using the prescription for protostellar jets presented in `Ostriker et al. (2001) <https://ui.adsabs.harvard.edu/abs/2001ApJ...557..443O/abstract>`_ and `Tabone et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018A%26A...614A.119T/abstract>`_. In this scenario, velocity variations within the beam of a highly collimated and highly supersonic jet induces de formation of internal working surface, from which jet material is ejected side ways, interacting with the ambient material. A bowshock shells produced from this interaction. Bowshockpy assumes full mixing between the jet and ambient material, that mass and momentum are conserved (negligible pressure gradients), and a negligible size of the internal working suface.

Workflow
--------

.. figure:: scheme_bowshockpy_text.png


    Workflow diagram of bowshockpy for computation of the bowshock model spectral cubes.


Below, we present a brief description of the workflow of ``bowshockpy``. The name of the main variables, as defined in :doc:`input parameters<inputparams>`, appear in parenthesis:

* From the mass and momentum conservation equations `Ostriker et al. (2001) <https://ui.adsabs.harvard.edu/abs/2001ApJ...557..443O/abstract>`_ and `Tabone et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018A%26A...614A.119T/abstract>`_, the morphology and kinematics of the bowshock shell can be obtained as a function of a few model parameters. These model parameters are characteristic length scale (L0), the distance between the working surface and the source (zj), the velocity of the internal working surfce (vj), velocity of the material surrounding jet (va), and the velocity at which the material is ejected from the internal working surface (v0).

* Given the total mass of the bowshock shell (mass), one can obtain its surface density. At this stage, we have all the parameters that define the model in its own reference frame. The rest of the workflow depends on the observer reference frame.

* In order to perform the mock observations, some parameters of dependent of the observer reference frame are used, mainly: the inclination angle of the bowshock axis (i), the observer distance to the source (distpc), the systemic velocity (vsys), and the position angle of the bowshock axis (PA). Together with some parameters defining the properties of the spectral cube, as the pixel size and channel width, ``bowshockpy`` computes, in projection, the mass of the bowshock shell at each pixel and velocity channel of the spectral cube. 

* ``bowshockpy`` can also calculate the intensity spectral cubes of a low-J rotational CO transtion. Assuming a CO abundance (XCO), ``bowshockpy`` calculates first the column densities at each pixel and channel of the cube. Given the excitation temperature (Tex) and assuming Local Thermodinamic Equilibrium, the opacities are computed. Finally, ``bowshockpy`` performs the radiative transfer in order to compute the intensities. 