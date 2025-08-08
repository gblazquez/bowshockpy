---
title: '`Bowshockpy`: A generator of synthetic observations for jet-driven bowshocks'
tags:
  - Python
  - astronomy
  - astrophysical jets
  - interstellar medium

authors:
  - name: Guillermo Blázquez-Calero
    orcid: 0000-0002-7078-2373
    corresponding: true
    affiliation: 1 
  - name: et al. 
    orcid: 0000-0000-0000-0000
    corresponding: true
    affiliation: 1 
affiliations:
 - name: Instituto de Astrofísica de Andalucía, CSIC, Glorieta de la Astronomía s/n, E-18008 Granada, Spain
   index: 1
date: 30 August 2025
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

<!--
Authors: Guillem Anglada, Mayra Osorio, Alex Raga, Sylvie Cabrit (?), Ruben Fedriani (?), Itziar de Gregorio, Jose Francisco Gomez, Gary Fuller, Noah Otten, Josep Maria, Rodrigo, Florin (?), Pablo (?), Oier (?)
-->

# Summary

Collimated jets are a common phenomenon in the Universe, appearing in a variety of astrophysical objects with accretion disks, such as black holes, planetary nebulae, and young stellar objects (YSOs). Astrophysical jets often exhibit knotty structures, which have long been suggested to be internal shocks resulting from velocity variations within the jet flow [@rees1978; @raga1992; @kofman1992]. In the case of YSO jets, interpreting the knots as "internal working surfaces" produced by a time-variable ejection provides an explanation for the multiple bow-shaped shocks typically observed in optical jets: the overpressure in the internal working surface drives material sideways, which then interacts with the ambient to form the bowshocks shells [@raga1990; @stone1993; @masson1993]. Current observations probing the molecular components of YSO jets --with high angular and spectral resolution-- <!-- with radio interferometers as the Atacama Large Milimmeter Telescope, are probing the molecular component of jets--> are now resolving shells connecting the jet knots, revealing their morphology and kinematics in great detail [@plunkett2015; @lee2017; @lopez-vazquez2024]. Modeling bowshock shells offers valuable insight into how jets interact with the environment [@blazquez-calero].

`Bowshockpy` is an open-source Python package that generates synthetic spectral cubes, position-velocity diagrams, and moment images of an analytical jet-driven bowshock model, using the prescription for YSO jets presented in @ostriker2001 and @tabone2018.  The software computes the intensities of low-J rotational transitions of the CO molecule, providing mock observations of the CO emission that radio telescopes as the Atacama Large Milimeter Array (ALMA) are able to detect at millimeter wavelengths.

<!--
TODO: Is the program generalizable for other CO rotational transition apart from CO(3-2)
-->


# Statement of need

Jets from YSOs are <!-- not a mere by-product of the star formation process, but are--> thought to play a crucial role in the formation of a star by removing angular momentum from the disk and, therefore, allowing accretion onto the forming star. Additionally, jets are invoked in order to explain the low star formation efficiency at both core and cloud scales [@frank2014]. Thus, characterizing the physical properties of YSO jets and their interaction with their surrounding medium is of major importance in order to understand the star formation process.  <!-- There are, however, some important unkowns; e.g., the launching mechanism of jets is still debated (resolving the launching zone, <1 au, is still not possible with the available instrumentation),  and the jet has been sometimes interpreted to be densest axial part of a radially extended wind [@wang2019] instead of being a truly narrow jet [@tafalla2017]. --> Recent observations at mm wavelengths with radio interferometers, specially ALMA, enable us to study in great detail the molecular component in jets. <!--which can shed light to these unkonws--><!--, nearest to the YSO ($\lesssim 5000$ au), and characterize these knots. At these scales,--> Sensitive observations at high and spectral resolution with ALMA reveal the presence of shell-like structures connecting the fast knots in the jet. While knots in the jets have been interpreted as internal working surfaces that expel jet material sideways by pressure forces [@tafalla2017], the shells could be bowshocks that arise from the interaction of this jet material with the surrounding gas [@plunkett2015; @lee2017; @lopez-vazquez2024]. Characterizing these bowshocks through models and numerical simulations [@ostriker2001; @tabone2018; @rabenanahary2022] can give information on the interaction of the jet with the surrounding gas, such as quantifying the velocity and mass-rate at which jet material is expeled sideways from the internal working surface, the mass-rate of ambient material incorporated into the bowshock and the ambient density [@blazquez-calero]. <!--, and constrain the launching mechanisms.--> <!-- provide valueable information of the jet properties, its surrounding ambient and their interaction.-->  

The necessity of `Bowshockpy` relies on the importance of providing a public open-source tool for modelling bowshock shells, with particular focus on the CO emission observable at millimeter wavelengths by radio telescopes as ALMA. Although there are numerous publications that have modeled many different aspects of bowshocks from YSO jets [@lee2000; @smith2003; @schultz2005; @correia2009; @burns2016; @rabenanahary2022; @rivera-ortiz2023], there is no public software designed for these puposes.  <!--that focus on modeling the molecular emission from bowshock shells.--> The analytical momentum-conserving bowshock jet-driven model implemented in `Bowshockpy` [@ostriker2001;@tabone2018], enables fast computation and, despite its assumptions (mainly, negligible pressure gradients and full mixing between jet and ambient material), its validity has been tested by more realistic hydrodynamcal simulations [@tabone2018].


# Code description

`Bowshockpy` computes synthetic spectral cubes, position-velocity diagrams, and moment images for an analytical jet-driven bowshock model, using the prescription for YSO jets presented in @ostriker2001 and @tabone2018. In this scenario, velocity variations within the beam of a highly collimated and highly supersonic jet induces the formation of internal working surfaces, from which jet material is ejected sideways. This jet material interacts with the ambient material, forming the bowshock shell. The analytical prescription implemented in `bowshockpy` assumes mass and momentum conservation (negligible preasure gradients), full mixing between the jet and ambient material, and a negligible size of the internal working suface.

These are the steps followed by `Bowshockpy` in order to obtain the CO intensity spectral cubes: 

- From the mass and momentum conservation equations, the morphology and kinematics of the bowshock shell can be obtained as a function of a few free parameters (see @ostriker2001 and @tabone2018). These model parameters are characteristic length scale<!--(L0)-->, the distance between the working surface and the source<!--(zj)-->, the velocity of the internal working surfce<!--(vj)-->, velocity of the material surrounding jet<!--(va)-->, and the velocity at which the material is ejected from the internal working surface <!--(v0)-->.

- Given the total mass of the bowshock shell, one can obtain its surface density. At this stage, we have all the parameters that define the model in its own reference frame. The rest of the workflow depends on the observer reference frame.

- In order to perform the mock observations, some parameters dependent of the observer reference frame are used, mainly: the inclination angle of the bowshock axis with respect to the line-of-sight<!--(i)-->, the observer distance to the source<!--(distpc)-->, the systemic velocity<!--(vsys)-->, and the position angle of the bowshock axis<!--(PA)-->. Together with some parameters defining the properties of the spectral cube, as the pixel size and channel width, bowshockpy computes, in projection, the mass of the bowshock shell at each pixel and velocity channel of the spectral cube.

- `Bowshockpy` can also calculate the intensities of low-J rotational CO transtions. Assuming a CO abundance<!--(XCO)-->, bowshockpy calculates first the column densities at each pixel and channel of the spectral cube. Given the excitation temperature<!--(Tex)--> and assuming Local Thermodinamic Equilibrium, the opacities are computed. Finally, `bowshockpy` performs the radiative transfer in order to compute the intensities.

The output spectral cubes, position velocity diagrams, and moments maps (integrated intensity, peak intensity, mean velocity field and vthe velocity dispersion) are saved in `FITS` format. 

The code is available at [GitHub](https://github.com/gblazquez/bowshockpy) and licensed under the MIT License. It can be install via PyPI, and it depends on the Python open-source libraries `NumPy` [@harris2020], `SciPy` [@virtanen2020], `Matplotlib` [@hunter2007], `Photutils` [@bradley2022], and `Astropy` [@astropy2013; @astropy2018; @astropy2022]. For more detailed information about `Bowshockpy` and its features, see the extensive documentation hosted at [ReadtheDocs](https://bowshockpy.readthedocs.io/), where examples and notebooks showing its usage are presented.

<!--
`Bowshockpy` is an implementation of the analytical momentum conserving bowshock model of @tabone2018, which has been validated with hydrodynamical numerical simulations. The software computes the observed column densities of the bowshock shell and, assuming local thermodynamic equilibrium conditions, simulates the CO intensities of low-J rotational transitions, providing mock observations of the CO emission that radio telescopes would be able to detect at milimeter wavelengths. This software provides spectral cubes, PV diagrams, and moment images of the bowshock shell fairly quick in an personal computer, allowing the modelization of observed bowshocks using just a few user-defined model and observer/instrumental parameters. To our knowledge, there is no software publicly available that allows the creation of synthetic observations of an analytical jet-driven bowshock model. We note that, although this software has been initially aimed to be used in jet-driven bowshocks in YSO objects, its applicability would be equally valid for molecular jets from protoplanetary nebulae.
-->

<!--
# Description

We summarize here the key principles and characteristics of the analytic, momentum-conserving bowshock model presented in Ref. @tabone2018, which we use as a basis for comparison with our data. Originally developed to describe the leading bowshock at the jet head[@masson1993; @ostriker2001], this model was recently extended by @tabone2018 to describe bowshocks formed by internal working surfaces (IWS) within a jet propagating into a slower-moving ambient medium.

In the framework of momentum conserving bowshock, velocity variations within a highly supersonic jet induces the formation of a two-shock structure called internal working surface [@raga1990]. Then, the overpressured shocked jet material is driven sideways, interacting with slower surrounding material, forming a curved bowshock. By modeling the bowshock as a stationary, thin shell of well-mixed material, its shape and velocity field can be derived self-consistently from the conservation of mass and momentum.


-->

# Acknowledgements

G.A., G.B.-C., IdG-M., G.A.F., J.F.G., M.O., acknowledge financial support from grants PID2020-114461GB-I00 and CEX2021-001131-S, funded by MCIN/AEI/10.13039/501100011033. G.B.-C., G.A.F., and M.O. acknowlege financial support from Junta de Andalucía (Spain) grant P20-00880 (FEDER, EU). G.B-C acknowledges support from grant PRE2018-086111, funded by MCIN/AEI/ 10.13039/501100011033 and by `ESF Investing in your future', and thanks ESO Science Support Discretionary Fund for their finantial support under the 2024 SSDF 06 project.  

# References
