# neutrino-winds
Neutrino-driven winds project
In order to successfully complete this assignment you need to turn in a project proposal to D2L on or before **11:59pm on Friday September 6**.

# <center> Neutrino-driven Wind Simulations </center>

<center>By Brian Nevins</center>

<img src="https://thumbor.forbes.com/thumbor/960x0/https%3A%2F%2Fblogs-images.forbes.com%2Fstartswithabang%2Ffiles%2F2017%2F05%2F1-T7LLvLFkGPm_iBOSl7af_g.jpg" width="60%">
<p style="text-align: right;">Image from: https://thumbor.forbes.com/thumbor/960x0/https%3A%2F%2Fblogs-images.forbes.com%2Fstartswithabang%2Ffiles%2F2017%2F05%2F1-T7LLvLFkGPm_iBOSl7af_g.jpg</p>

---
# Overview

This project will be focused on the subfield of nuclear astrophysics that studies proto-neutron stars. Proto-neutron stars are the core that is left behind when a massive star goes supernova. Proto-neutron stars undergo a brief period of stabilization and cooling, during which some of its mass is carried off into space by winds driven by the release of large quantities of neutrinos in the core. This matter is very neutron-rich, which allows for many interesting nuclear reactions and element synthesis to occur. The interior of the proto-neutron star is very active, and the convection inside creates gravitational waves which can also contribute to the nuclear reactions that take place in this wind. 

The behavior of these winds determines the kinds of nuclei that are formed, and this behavior is governed by a set of coupled differential equations. Computation is used to explore the behavior of these equations given certain initial conditions. The objective of this project is to determine which sets of boundary conditions result in a 'trans-sonic' wind - a wind faster than the speed of sound in ejected stellar media. Ultimately, we hope to use this information to determine the temperature and density of the neutrino-driven winds, which will allow us to determine the kinds of nuclei that are likely to be formed in them. We also hope to explore how the additional heating caused by gravitational waves inside the star affects the nuclear reactions that take place.

---
# Program Description

Some code exists that explores the steady-state behavior of the neutrino-driven winds in the isothermal approximation. I will be working to implement this in Python. This will involve implementing an efficient iterative method for evaluating the behavior of coupled differential equations, and solving these equations numerous times to determine the initial conditions that allow for the desired outcome. Ideally, I will also be able to investigate variations of the differential equations that include other heating effects. It will be helpful to generate visuals that effectively communicate the behavior and structure of these winds, and the boundary conditions required for useful data.

---
# Project Goals

Short-term:<br>
Implement a functioning iterative DE approximation method <br>
Generate useful visuals for DE behaviors

Mid-term:<br>
Efficiently implement a method for finding sets of initial conditions that produce the desired behavior in the DEs<br>
Improve accuracy of DE approximations without overly compromising efficiency

Long-term:<br>
Implement a minimization algorithm suggested by my advisor that may allow for better simulations of wind behavior

---
# Anticipating Challenges  

Understanding the physics of these winds well enough to implement the equations<br>
Learning to generate the visuals I need, especially in real-time<br>
Running large numbers of DEs accurately without taking too much time