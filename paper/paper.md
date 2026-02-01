---
title: '`BowshockPy`: A generator of synthetic observations for jet-driven bowshocks'
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
  - name: Guillem Anglada
    orcid: 0000-0002-7506-5429
    corresponding: false
    affiliation: 1 
  - name: Alejandro C. Raga
    orcid: 0000-0002-0835-1126
    corresponding: false
    affiliation: 2 
  - name: Sylvie Cabrit
    orcid: 0000-0002-1593-3693
    corresponding: false
    affiliation: "3, 4"
  - name: Mayra Osorio
    orcid: 0000-0002-6737-5267
    corresponding: false
    affiliation: 1 
  - name: José F. Gómez 
    orcid: 0000-0002-7065-542X
    corresponding: false
    affiliation: 1 
  - name: Gary A. Fuller
    orcid: 0000-0001-8509-1818
    corresponding: false
    affiliation: 5 
  - name: Ruben Fedriani
    orcid: 0000-0003-4040-4934
    corresponding: false
    affiliation: 1 
  - name: Itziar de Gregorio-Monsalvo
    orcid: 0000-0003-4518-407X
    corresponding: false
    affiliation: 7 
  - name: Noah Otten
    orcid: 0000-0002-2530-4137
    corresponding: false 
    affiliation: 8 
  - name: Florin Placinta-Mitrea
    orcid: 0009-0009-7062-3629
    corresponding: false
    affiliation: 1 
  - name: Pablo Santo-Tomás
    orcid:  0009-0003-5228-2804
    corresponding: false
    affiliation: 1 
  - name: Oier Baraibar
    orcid: 0009-0003-1960-2327
    corresponding: false
    affiliation: 1 
  - name: Josep M. Masqué
    orcid: 0000-0001-9553-6614
    corresponding: false
    affiliation: 9 
  - name: Rodrigo D. Garduño
    orcid: 0009-0006-9272-308X
    corresponding: false
    affiliation: 10

affiliations:
 - name: Instituto de Astrofísica de Andalucía, CSIC, Glorieta de la Astronomía s/n, E-18008 Granada, Spain
   index: 1
 - name: Instituto de Ciencias Nucleares, Universidad Nacional Autónoma de México, Apartado Postal 70-543, 04510 Ciudad de México, Mexico
   index: 2
 - name: Observatoire de Paris --- PSL University, Sorbonne Université, CNRS, 75014 Paris, France
   index: 3
 - name: Univ. Grenoble Alpes, CNRS, IPAG, 38000 Grenoble, France
   index: 4
 - name: Jodrell Bank Centre for Astrophysics, Department of Physics and Astronomy, The University of Manchester, Oxford Road, Manchester M13 9PL, UK
   index: 5
 - name: I. Physikalisches Institut, University of Cologne, Zülpicher Str. 77, 50937 Küoln, Germany
   index: 6
 - name: European Southern Observatory, Alonso de Córdova 3107, Casilla 19, Vitacura, Santiago, Chile
   index: 7
 - name: Department of Physics, Maynooth University, Maynooth, Co. Kildare, Ireland
   index: 8
 - name: Institut de Ciències del Cosmos (ICCUB), Universitat de Barcelona (UB), Martí de Franquès 1, E-08028 Barcelona, Spain
   index: 9
 - name: Departamento de Astronomía, Universidad de Guanajuato, Apartado Postal 144, 36000 Guanajuato, México
   index: 10


date: 30 August 2025
bibliography: paper.bib

---


# Summary

Collimated jets are a common phenomenon in the Universe, appearing in a variety of astrophysical objects with accretion disks, such as black holes, planetary nebulae, post-AGB stars, and young stellar objects (YSOs). Astrophysical jets often exhibit knotty structures, which have long been suggested to be internal shocks resulting from velocity variations within the jet flow [@rees1978; @raga1992; @kofman1992]. In fact, numerical simulations show that supersonic variations in the ejection velocity lead to the formation of two-shock structures known as "internal working surfaces", which travel downstream of the jet flow [@biro1994]. The overpressure in the internal working surface drives material sideways that interacts with the ambient, producing the bow-shaped shocks typically observed in jets from YSO [@raga1990; @stone1993; @masson1993]. Modeling bowshock shells offers valuable insight into how jets interact with the ambient medium, enabling us to characterize jet properties that are essential for understanding how stars form [@blazquezcalero2025].

`BowshockPy` is an open-source Python package that generates synthetic spectral cubes, position-velocity diagrams, and moment images of an analytical jet-driven bowshock model, using the prescription for YSO jets presented in @ostriker2001 and @tabone2018. `BowshockPy` assumes a thin bowshock shell whose morphology and kinematics are determined by the mass and momentum conservation (ignoring pressure gradients), considering full mixing between the jet and ambient material as well as a negligible size of the internal working surface. The software calculates the column densities of the bowshock shell and can determine the intensities of low-J rotational transitions of linear molecules (such as CO). The intensities are obtained using the rigid rotor approximation (valid for low-J transitions, where vibrational excitation and centrifugal distortion effects are negligible), and assuming local thermodynamic equilibrium. This provides mock observations of the molecular line emission that radio telescopes as the Atacama Large Millimeter Array (ALMA) are able to image at millimeter wavelengths.


# Statement of need

Jets from YSOs are thought to play a crucial role in the formation of a star by removing angular momentum from the disk and, therefore, allowing accretion onto the forming star. Additionally, jets are invoked in order to explain the low star formation efficiency at both core and cloud scales [@frank2014]. Thus, characterizing the physical properties of YSO jets and their interaction with their surrounding medium is of major importance in order to understand the star formation process. Recent observations at mm wavelengths with radio interferometers, specially ALMA, enable us to study in great detail the molecular component in jets. Sensitive observations at high spatial and spectral resolution with ALMA reveal the presence of shell-like structures connecting the fast knots in the jet. While knots in the jets have been interpreted as internal working surfaces that eject jet material sideways by pressure forces [@santiago-garcia2009], the shells could be bowshocks that arise from the interaction of this jet material with the surrounding gas [@plunkett2015; @lee2017; @lopez-vazquez2024]. Characterizing these bowshocks through models and numerical simulations [@ostriker2001; @tafalla2017; @tabone2018; @rabenanahary2022; @rivera-ortiz2023; @tafalla2026] can provide quantitative information on the interaction of the jet with the surrounding gas, such as the velocity and mass-rate at which jet material is ejected sideways from the internal working surface, the mass-rate of ambient material incorporated into the bowshock, and the ambient density [@blazquezcalero2025]. This information helps to understand the dynamical properties of jets and how it injects mass and momentum to the environment, both of which are crucial for understanding how stars forms.

The necessity of `BowshockPy` relies on the importance of providing a public open-source tool for modeling bowshock shells, with particular focus on the line emission of low-J rotational transitions of linear molecules (such as CO), observable at millimeter wavelengths by radio interferometers as ALMA.


# State of the field

Although there are numerous publications that have modeled many different aspects of bowshocks from YSO jets [@lee2000; @smith2003; @schultz2005; @correia2009; @burns2016; @rabenanahary2022; @rivera-ortiz2023; @tafalla2026], there is no public software designed for these purposes. The simple analytical momentum-conserving bowshock model implemented in `BowshockPy` [@ostriker2001;@tabone2018], which assumes a thin shell of fully mixed jet and ambient material, enables fast computation. The model results have been compared with more detailed hydrodynamical simulations, obtaining a good agreement in the bowshock morphology and kinematics [@tabone2018]. 


# Software design

`BowshockPy` is a Python 3 software that uses a modular architecture, designed with the purpose of generating in a simple and quick way model images of bowshocks that can be directly compared with observations. `BowshockPy` computes synthetic spectral cubes, position-velocity diagrams, and moment maps (integrated intensity, peak intensity, mean velocity field, and velocity dispersion) for an analytical jet-driven bowshock model, based on the prescription for YSO jets presented in @ostriker2001 and @tabone2018. These are the steps followed by `BowshockPy` in order to obtain the line intensity spectral cubes:

- From the mass and momentum conservation equations, the morphology and kinematics of the bowshock shell can be obtained as a function of a few free parameters [@ostriker2001;@tabone2018]. These model parameters are a characteristic length scale, the distance between the working surface and the source, the velocity of the internal working surface, the velocity of ambient material surrounding the jet, and the velocity at which the material is ejected from the internal working surface. The surface density at each point of the bowshock is computed as a function of the shell integrated mass [@blazquezcalero2025]. At this stage, we have all the parameters that define the model in its own reference frame. The rest of the workflow depends on the observer reference frame.

- In order to perform the mock observations, some parameters dependent of the observer reference frame are used, mainly: the inclination angle of the bowshock axis with respect to the line-of-sight, the observer distance to the source, the systemic velocity of the ambient cloud, and the position angle of the bowshock axis.

- Together with some parameters defining the properties of the spectral cube, such as the pixel size and channel width, `BowshockPy` computes the mass of the bowshock shell at each pixel and velocity channel of the spectral cube.

- From the masses at each pixel and channel of the spectral cube, `BowshockPy` computes the column densities of the emitting molecule assuming a given abundance. In addition, it can calculate the opacities under local thermodynamic equilibrium for a low-J rotational transition of a linear molecule in the rigid rotor approximation (that is, under the assumption of negligible vibrational excited states and negligible centrifugal distortion effects), and perform the radiative transfer to obtain the intensities. If the user needs a different model to determine the intensities, `BowshockPy` allows them to apply custom models to the calculated column densities.

There are two ways to utilize `BowshockPy`: It can either be run from the terminal using an input file containing all the model parameters, or it can be imported as a package to use its functions and classes. The computed spectral cubes, position velocity diagrams, and moments maps are saved in `FITS` format. For more detailed information about `BowshockPy` and its features, see the extensive documentation hosted at [ReadtheDocs](https://bowshockpy.readthedocs.io/), where examples and a notebook tutorial showing its usage are presented.

The code is available at [GitHub](https://github.com/gblazquez/bowshockpy) and licensed under the MIT License. It can be installed via PyPI, and it depends on the Python open-source libraries `NumPy` [@harris2020], `SciPy` [@virtanen2020], `Matplotlib` [@hunter2007], `Photutils` [@bradley2022], and `Astropy` [@astropy2013; @astropy2018; @astropy2022]. 


# Research impact Statement

`BowshockPy` has proven to be highly valuable for modeling jet-driven bowshock shells. It has been used in @blazquezcalero2025, which showed that the bowshock model is able to reproduce the morphology and kinematics of the observed shell-like structures in a molecular jet. Its capability to generate images that can be directly compared with observations at a low computational cost makes it a very useful tool for the scientific community specialized in jets from YSO, particularly for interpreting sensitive radio interferometric observations with high spatial and spectral resolution. 


# AI usage disclosure

No AI tools were used in the creation of this software. Generative AI tools were used only occasionally to correct the spelling and improve the clarity of some parts of the text in this paper.


# Acknowledgments

G.A., G.B.-C., IdG-M., G.A.F., J.F.G., M.O., R.F., F.P.M., O.B., P.S.-T., acknowledge financial support from grants PID2023-146295NB-I00PID and CEX2021-001131-S, funded by MCIN/AEI/10.13039/501100011033. G.B-C acknowledges support from grant PRE2018-086111, funded by MCIN/AEI/10.13039/501100011033 and by `ESF Investing in your future', and thanks ESO Science Support Discretionary Fund for their finantial support under the 2024 SSDF 06 project. S.C. acknowledges support by the Thematic Action “Physique et Chimie du Milieu Interstellaire” (PCMI) of INSU Programme National “Astro”, with contributions from CNRS Physique & CNRS Chimie, CEA, and CNES. N. Otten acknowledges funding from the Maynooth University Graduate Teaching Scholarship and Taighde Éireann (Research Ireland) under the RI-ESO Studentship Agreement for this work. F.P.M. acknowledges support from grant PRE2021-100926, funded by MCIN/AEI/10.13039/501100011033. P.S.-T. acknowledges financial support from project ref. AST22_00001_Subp X with funding from the European Union - NextGenerationEU; the Spanish Ministry of Science, Innovation and Universities; the Spanish Recovery, Transformation and Resilience Plan; the Department of University, Research and Innovation of the Andalusian Regional Government and Consejo Superior de Investigaciones Científicas. O.B. acknowledges support from grant PREP2023-001438, funded by MCIN/AEI/10.13039/501100011033.  

# References
