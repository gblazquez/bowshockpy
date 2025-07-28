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
Authors: Guillem Anglada, Mayra Osorio, Sylvie Cabrit (?), Benoit Tabone (?), Ruben Fedriani (?), Alejandro López-Vazquez (?), Itziar de Gregorio, Jose Francisco Gomez, Gary Fuller, Noah Otten, Josep Maria, Rodrigo, Florin (?), Pablo (?), Oier (?)
-->

# Summary
<!--
# Introduction
# Background
-->
<!--
Possibility 1:

Ejections in star formation. Ejection mechanisim is unkown. Molecular jets. X-wind vs bowshock. Bowshock faint. Properties of the interaction can be obtain.

Possibility 2:

Collimated jets in star formation. Optical/IR are low resolution. Molecular component in radio: high resolution, enabling comparison with models. Apart from the molecular jets, bowshock wings as nested shells are being observed with radio interferometers.
-->

Collimated jets are a common phenomenon in the Universe, appearing in a variety of astrophysical objects with accretion disks, such as black holes, planetary nebulae, and young stellar objects (YSOs). These astrophysical jets sometimes present knotty structures, which have long been suggested to be the result of a velocity variations of the flow within the jet beam [@rees1978]. In the case of jets from YSOs, interpreting the knots as "internal working surfaces" produced by time-variable ejection, provided an explanation for the multiple bow-shaped shocks typically observed in Herbig-Haro jets: the overpressure in the internal working surface drives material sideways, which interacts with the ambient to form the bowshocks shells [@raga1990; @stone1993; @masson1993]. Current observations probing the molecular components of jets <!-- with radio interferometers as the Atacama Large Milimmeter Telescope, are probing the molecular component of jets--> from YSOs are detecting bowshock shells with sufficient high angular and spectral resolution to reveal in details of their morphology and kinematics. The characterization of these bowshock using models provides valueable informaton of the interaction of the jet with the environment (Blázquez-Calero et al., under review). `Bowshockpy` is a Python package that generates synthetic spectral cubes, PV diagrams, and moment images for a simple analytical jet-driven bowshock model, using the prescription for protostellar jets presented in @tabone2018.

<!--
TODO: Is the program generalizable for other CO rotational transition apart from CO(3-2)
-->

<!--
Along with observer parameters as the inclination angle, the projected morphology and kinematics are obtained. The surface density is computed, as well as the intensity .
-->

# Statement of need
<!--
time-dependent ejections 
variations in the flow velocity/ ejection velocty within the jet beam
-->
<!--
In the case of YSO, Hypersonic, collimated protostelar mass loss appearas to be a ubiquitous aspect of the star formaton process.
The suggestion that the knotty structures in astro-
physical jets could be the result of a time-dependent
ejection was first made in the context of extragalac-
tic jets (see, e.g., Rees 1978; Wilson 1984; Roberts
1986). However, the theory of variable jets has
been mostly developed and applied in the context
of Herbig-Haro (HH) jets from young stars.
Raga et al. (1990) apparently first pointed out
in an explicit way that the structures observed in
HH jets could be easily modeled as “internal working
surfaces” produced by an ejection velocity variabil-
ity with a hypersonic amplitude (though the general
idea that HH knots are the result of a variability of
the ejection hovers around in the literature of the late
1980’s)
--> 
<!--
the star formation process is accompanied with the ejection of matter in the form of highly supersonic jets [@frank2014].

Jets from YSOs are thought to play an important role in the formation of a star by removing angular momentum from disk and, therefore, allowing accretion onto the central young star.
-->

<!-- say that YSO jets are supersonic and radiative?-->

<!--
Jets from YSOs play an important role in the formation of a star by removing angular momentum from disk and, therefore, allowing accretion onto the central young star[@frank2014]. However, the launching mechanism of jets is still debated, since resolving the launching zone (<1 au) is still not possible with the available instrumentation. 


Thus, indirect ways of constraining the ejection mechanism are needed. One way to constrain the launching properties is the characterization of the jets properties through its interaction with the surrounding gas, which can give insight into the mass-loss rate.

The interaction of jet internal working surface with the surrounding medium can give as insight of the mass-loss rate. Thus, the characterization of jets through its interaction with the environment is very important. 
-->

Jets from YSOs are <!-- not a mere by-product of the star formation process, but are--> thought to play a crucial role in the formation of a star by removing angular momentum from the disk and, therefore, allowing accretion onto the forming star. Additionally, jets are invoked in order to explain the low star formation efficiency at both core and cloud scales [@frank2014]. Thus, characterizing the physical properties of protostellar jets and their interaction with their surrounding medium is of major importance in order to understand the star formation process.  <!-- There are, however, some important unkowns; e.g., the launching mechanism of jets is still debated (resolving the launching zone, <1 au, is still not possible with the available instrumentation),  and the jet has been sometimes interpreted to be densest axial part of a radially extended wind [@wang2019] instead of being a truly narrow jet [@tafalla2017].  --> Recent observations at mm wavelengths with radio interferometers, specially the Atacama Large Milimeter Array (ALMA), enable us to study in great detail the molecular component in jets (mainly in CO and SiO emisson). <!--which can shed light to these unkonws--><!--, nearest to the YSO ($\lesssim 5000$ au), and characterize these knots. At these scales,--> Sensitive observations at high and spectral resoltion with ALMA reveal the presence of shell-like structures connecting the fast knots in the jet. While knots in the jets are usually interpreted as internal working surfaces [@tafalla2017] that expel jet material sideways by pressure forces, the shells are thought to be bowshocks that arise from the interaction of this jet material with the surrounding gas [@plunkett2015; @lee2017]. Characterizing these bowshocks through models an numerical simulations [@ostriker2001; @tabone2018] can give information on the interaction of the jet with the surrounding gas, such as quantifying the velocity and mass-rate at which jet material is expeled sideways from the internal working surface, the mass-rate of ambient material incorporated into the bowshock and the ambient density (Blázquez-Calero et al., under review). <!--, and constrain the launching mechanisms.--> <!-- provide valueable information of the jet properties, its surrounding ambient and their interaction.-->  

<!-- 
`Bowshockpy` is a Python package that generates synthetic spectral cubes, PV diagrams, and moment images for a simple analytical jet-driven bowshock model, using the prescription presented in @tabone2018 for protostellar jets. The code is an implementation of @tabone2018 prescription for a protostellar jet propagating in a surrounding disk wind.  In this framework, velocity variations within the jet beam induces the formation of internal working surfaces, from which the jet material is ejected sideways. These jet material interacts with the surrounding medium, forming a momentum conserving bowshock shell of well mixed ambient and jet material. `Bowshockpy` computes the morphology, kinematics, and surface density of the bowshock shell using a few user-defined model parameters. Then, assuming some user-defined observer and instrumental properties, this software simulates the observed column densities and, under local thermodynamic equilibrium conditions, computes the CO intensities of low-J rotational transitions, providing mock observations of the CO emission that mm radio telescopes are able to detect.
--> 

<!--
, providing an explanation of the multiple bow shock structures observed in some jets from YSOs.

Momentum conserving bowshock models are found in literature [@ostriker2001; @tabone2018]. Nonetheless, it has not been until recently that, by the advent of mm radio interferometers as ALMA, we obtained observations with enough angular and spectral resolution, sensitive enough to detect and model bowshocks (Blázquez-Calero et al., under rev.), mainly within the <5000 au. When compared to observations, the characterization of bowshocks can give information on the interaction with the ambient medium / entrainment process / and can even elucidate the launching mechanism (tafalla vs wang). 

-->

<!--
Jets from YSOs play a key role in the formation of a star by removing angular momentum from the star/disk system, however its launching mechanism is still debated. Resolving the launching zone (<1 au) is still not possible with the available instrumentation, so indirect ways are needed. 

 Bowshock shells are 
-->

`Bowshockpy` is an implementation of the analytical momentum conserving bowshock model of @tabone2018, which has been validated with hydrodynamical numerical simulations. The software computes the observed column densities of the bowshock shell and, assuming local thermodynamic equilibrium conditions, simulates the CO intensities of low-J rotational transitions, providing mock observations of the CO emission that radio telescopes would be able to detect at milimeter wavelengths. This software provides spectral cubes, PV diagrams, and moment images of the bowshock shell fairly quick in an personal computer, allowing the modelization of observed bowshocks using just a few user-defined model and observer/instrumental parameters. <!--, and its applicability is two-fold. First, enables the modelization of observed bowshocks using a just few user-defined model and observer/instrumental parameters.  Second, it can be used as a first approach for tailored computational expensive numerical magneto-hydrodynamical simuations. --> To our knowledge, there is no software publicly available that allows the creation of synthetic observations of an analytical jet-driven bowshock model. We note that, although this software has been initially aimed to be used in jet-driven bowshocks in YSO objects, its applicability would be equally valid for molecular jets from protoplanetary nebulae.

<!--
Analytical model that is computed quickly. Visualize and quickly compare with observations. This software is of scientific important since it enables to characterize the interaction between jets and environment. Also, guess parameters for time consuming MHD simulations.

. Moreover, the modelization of
bowshocks could potentially distinguiwish between ejection mechanisms
(tafalla2017, wang2019).

and create synthetic observations. 
-->

<!--
# Description

We summarize here the key principles and characteristics of the analytic, momentum-conserving bowshock model presented in Ref. @tabone2018, which we use as a basis for comparison with our data. Originally developed to describe the leading bowshock at the jet head[@masson1993; @ostriker2001], this model was recently extended by @tabone2018 to describe bowshocks formed by internal working surfaces (IWS) within a jet propagating into a slower-moving ambient medium.

In the framework of momentum conserving bowshock, velocity variations within a highly supersonic jet induces the formation of a two-shock structure called internal working surface [@raga1990]. Then, the overpressured shocked jet material is driven sideways, interacting with slower surrounding material, forming a curved bowshock. By modeling the bowshock as a stationary, thin shell of well-mixed material, its shape and velocity field can be derived self-consistently from the conservation of mass and momentum.


Summary: Few parameters that define the bowshock, we obtain the CO spectral
cube, pv's, and moments.

- Explain the bowshock model. Foundations (references). Morphology and
  kinematics given in Ostriker and Tabone. Parameters that define a bowshock. 
- Surface density (ref of your paper?)
- Mass in each pixel. CIC interpolation
- Once we have the mass in each cell of the spectral cube through equation, we can calculate the intensity of the line of interest, assuming the excitation properties and performing the radiative transfer. In this thesis, we are interested in the CO emission from a bowshock model, assuming LTE conditions and perform for the radiative transfer (eqs from ). In order to compare it with radio observations, we convolved the model images with the synthetized beam of the observations.
- Outputs: Cube, pv's and moments, but also important parameters such as mdot0,
  mdotamb, and the ambient density.
-->

<!--
that results from the mass and $(x^*,r)$-momentum conservation equations:
\begin{eqnarray}
	{\dot m} & = & {\dot m}_0+\pi r_b^2 \rho_{\rm amb}(v_{\rm jet}-v_{\rm amb})=2\pi r_b \sigma v_t\,, \label{eq:mcon} \\ 
	{\dot \Pi}_{x^*} & = & \pi r_b^2\rho_{\rm amb}(v_{\rm jet}-v_{\rm amb})^2={\dot m} v_{x^*}\,, \label{eq:xcon} \\
  {\dot \Pi}_r & = & {\dot m_0}v_0={\dot m}v_r\,,
  \label{eq:rcon}
\end{eqnarray}
where ${\dot m}$, ${\dot \Pi}_{x^*}$ and ${\dot \Pi}_r$ are the mass, $x^*$-momentum and $r$-momentum rates flowing along the thin shell up to a given value of $x^*$, and $v_t$, $v_{x^*}$ and $v_r$ are the components of the velocity of the well mixed material within the shell along the shell surface, and along the $x^*$- and $r$-axes, respectively. Finally, $\sigma$ (see the last term of \autoref{eq:mcon} is the surface density of the thin shell.

\begin{eqnarray}
	v_{x^*} & = &\frac{\pi r_b^2\rho_{\rm amb}(v_{\rm jet}-v_{\rm amb})^2}{\dot m_0+\pi\rho_{\rm amb}(v_{\rm jet}-v_{\rm amb})r_b^2}\,, \label{vx} \\ 
	v_r & = & \frac{{\dot m}_0v_0}{\dot m_0+\pi\rho_{\rm amb}(v_{\rm jet}-v_{\rm amb})r_b^2}\,.
  \label{vr}
\end{eqnarray}
integrated to obtain the shape of the bowshock wings:
\begin{equation}
	r_b(x^*)=\left(L_0^2\, x^*\right)^{1/3},
  \label{rb}
\end{equation}
with
\begin{equation}
	L_0\equiv \sqrt{\frac{3{\dot m}_0 v_0}{\pi\rho_{\rm amb}(v_{\rm jet}-v_{\rm amb})^2}}\,.
  \label{l0}
\end{equation}

Finally, we consider that the velocity along the thin shell can be written as $v_t=v_{x^*}\cos\alpha+v_r\sin\alpha$, and use equations (\ref{eq:mcon}) and (\ref{eq:vx})-(\ref{drb}) to calculate
the surface density of the shell as
\begin{equation}
	\sigma=\frac{1}{2}~\rho_{\rm amb}\cos\alpha\left(\gamma\tan\alpha + 1\right)^2 r_b\,,
	\label{sig}
\end{equation}
with


- The free parameters are...

- We provide a Documentation with examples of the code usage.

-->

<!--
# Figures

Figures can be included like this:

![Caption for example figure.\label{fig:example}](scheme_bowshockpy.pdf){ width=100% }

and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }
-->

# Acknowledgements

We acknowledge the significant contribution of Alex Raga to this work, who made a preliminary version of this software in Fortran. G.A., G.B.-C., IdG-M., A.K.D.-R., G.A.F., J.F.G., M.O., acknowledge financial support from grants PID2020-114461GB-I00 and CEX2021-001131-S, funded by MCIN/AEI/10.13039/501100011033. G.B.-C., G.A.F., and M.O. acknowlege financial support from Junta de Andalucía (Spain) grant P20-00880 (FEDER, EU). G.B-C acknowledges support from grant PRE2018-086111, funded by MCIN/AEI/ 10.13039/501100011033 and by `ESF Investing in your future', and thanks ESO Science Support Discretionary Fund for their finantial support under the 2024 SSDF 06 project.  

# References
