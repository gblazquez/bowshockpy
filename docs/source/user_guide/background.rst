Background
====================

The bowshock model
------------------

``BowshockPy`` computes synthetic spectral cubes, PV diagrams, and moment images for a simple analytical jet-driven bowshock model, using the prescription for protostellar jets presented in `Ostriker et al. (2001) <https://ui.adsabs.harvard.edu/abs/2001ApJ...557..443O/abstract>`_ and `Tabone et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018A%26A...614A.119T/abstract>`_. In this scenario, velocity variations within the beam of a highly collimated and highly supersonic jet induces the formation of internal working surfaces, from which jet material is ejected sideways. This jet material interacts with the ambient material, forming the bowshock shell. The analytical prescription implemented in ``BowshockPy`` assumes a thin bowshock shell whose morphology and kinematics are determined by the mass and momentum conservation (with negligible pressure gradients), full mixing between the jet and ambient material, and a negligible size of the internal working surface.

.. 
   Although the model was focused on bowshocks from protostellar jets, we note that it could also work for jets associated to proto-planetary nebulae.

Workflow
--------

.. figure:: images/scheme_bowshockpy_text.png


    Workflow diagram of ``BowshockPy`` for computation of the bowshock model spectral cubes.


Below, we present a brief description of the workflow of ``BowshockPy``. The name of the model parameters, as defined in :doc:`input parameters<inputparams>`, appear in parenthesis:

* From the mass and momentum conservation equations, the morphology and kinematics of the bowshock shell can be obtained as a function of a few free parameters (see `Ostriker et al. (2001) <https://ui.adsabs.harvard.edu/abs/2001ApJ...557..443O/abstract>`_ and  `Tabone et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018A%26A...614A.119T/abstract>`_). These model parameters are characteristic length scale (L0), the distance between the working surface and the source (zj), the velocity of the internal working surface (vj), velocity of the material surrounding jet (va), and the velocity at which the material is ejected from the internal working surface (v0). The surface density at each point of the bowshock is computed as a function of the shell integrated mass (mass) `Blazquez-Calero et al. (2025) <https://ui.adsabs.harvard.edu/abs/2025NatAs.tmp..254B/abstract>`_. At this stage, we have all the parameters that define the model in its own reference frame. The rest of the workflow depends on the observer reference frame.

* In order to perform the mock observations, some parameters dependent of the observer reference frame are used, mainly: the inclination angle of the bowshock axis with respect to the line-of-sight (i), the observer distance to the source (distpc), the systemic velocity (vsys), and the position angle of the bowshock axis (PA). 

* Together with some parameters defining the properties of the spectral cube, as the pixel size and channel width, ``BowshockPy`` computes, in projection, the mass of the bowshock shell at each pixel and velocity channel of the spectral cube. 

* ``BowshockPy`` then computes the column densities of the emitting molecule assuming an abundance relative to the molecular hydrogen (abund). In addition, it can calculate the opacities under local thermodynamic equilibrium for a low-J rotational transition of a linear molecule in the rigid rotor approximation (that is, under the assumption of negligible vibrational excited states and negligible centrifugal distortion effects), and perform the radiative transfer assuming an excitation temperature (Tex) to obtain the intensities. If the user requires a different modeling to obtain the intensities, ``BowshockPy`` allows them to apply a custom model to the computed column densities.

   
References
----------

- Tabone, B., Raga, A., Cabrit, S. & Pineau des Forêts, G. "Interaction between a pulsating jet and a surrounding disk wind. A hydrodynamical perspective." Astron. Astrophys. 614, A119 (2018).

- Ostriker, E. C., Lee, C.-F., Stone, J. M. & Mundy, L. G. "A Ballistic Bow Shock Model for Jet-driven Protostellar Outflow Shells". Astrophys. J. 557, 443–450 (2001).

- Blazquez-Calero, G., Anglada, G., Cabrit, S., Osorio, M., et al. "Bowshocks driven by the pole-on molecular jet of outbursting protostar SVS 13". Nature Astronomy (2025). doi:10.1038/s41550-025-02716-2

.. _Tabone et al. (2018): https://ui.adsabs.harvard.edu/abs/2018A%26A...614A.119T/abstract
.. _Ostriker et al. (2001): https://ui.adsabs.harvard.edu/abs/2001ApJ...557..443O/abstract
.. _Blazquez-Calero et al. (2025): https://ui.adsabs.harvard.edu/abs/2025NatAs.tmp..254B/abstract
 
