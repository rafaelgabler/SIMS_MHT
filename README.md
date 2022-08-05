# SIMS_MHT

Developed by Rafael Gabler Gontijo, PhD

- Numerical simulation of multibody suspensions of magnetic particles using Langevin and Stokesian Dynamics. 
- Optimized version for magnetic hyperthermia simulations.

This code solves the translational and rotational motion of N solid, spherical, magnetic particles, suspended in a viscous liquid and subjected to an external magnetic field. 

The code is composed by the following files:

entrada.dat = configuration file in which the user sets simulation parameters 

entrada.f90 = fortran module that reads the configuration file and stores its information in terms of logical, integer and real variables

sims.f90 = a small file that simply starts the simulation process and count the computational time

principal.f90 = main file which calls all the subroutines in a logical sequence in order to perform the simulation set in the configuration file

funcoes.f90 = fortran module containing all the subroutines used in the code

variables.f90 = fortran module which defines all the variables used in the code

saida.f90 = performs several statistical analysis based on velocities, orientations and position of the particles

How it works?

Sims was designed to perform simultaneous numerical experiments of N magnetic particles immersed in a viscous liquid during sedimentation or shear induced motion. The particles can be either neutrally buoyant or denser than the base fluid. They can undergo a steady state or either an unsteady oscillatory shear field. They can also be in the presence or in the absence of an applied magnetic field. When in the presence of an applied field, this field can be either steady-state or time dependent. In the case of time dependent field excitation the user can choose an oscillatory 1D field, a rotating 2D field, a non-linear 1D Duffing excitation or even a 1D oscillatory beating pattern excitation. The user can also sets other options, such as: compute or neglect Brownian effects, change the initial configuration of the particles to a set of ordered spheres or to a random dispersion of spheres inside a spherical aggregate-like structure, mix magnetic with non-magnetic particles and compute only velocity fluctuations. 

The user defines what scenario should be simulated by changing the data displayed in the configuration file 'entrada.dat'. There, the user defines the values of logical variables (true or false) by basically answering some questions regarding what is going to be simulated. Then, the user can set numerical data regarding the number of particles and the number of independent realizations, the volume fraction of particles, the aspect lenght of the simulation box and the number of physical and reciprocal lattice cells used to compute long range interactions, which are relevant only when we compute hydrodynamic interactions and dipolar magnetic torques in concentrated regimes where the volume fraction of particles is higher than 10%. The number 125 for this option should work fine. The user can also set the frequencies of the applied shear rate and magnetic field, when it is the case. 

Informations regarding the intensity of the applied field, magnitude of magnetic moments, size of the particles, viscosity of the surrounding liquid and local temperature are expressed in terms of non-dimensional physical parameters, namely:

St = Stokes number

Str = Rotational Stokes number

Pe = Péclet number

alpha = dimensionless magnetic field

lambda = dimensionless intensity of the magnitude of the magnetic moments of the particles

For more information regarding the mathematical formulation used in SIMS, examples and a better physical interpretation of the dimensionless parameters mentioned above, the user should feel free to consult the following references.

- GONTIJO, R.G.; CUNHA, F.R. . Dynamic numerical simulations of magnetically interacting suspensions in creeping flow. Powder Technology (Print), v. 279, p. 146-165, 2015.

- GONTIJO, R.G.; Malvar, S. ; CUNHA, F.R. . Magnetic particulate suspensions from the perspective of a dynamical system. Powder Technology (Print), v. 297, p. 165-182, 2016.

- GONTIJO, R.G.. A numerical perspective on the relation between particle rotational inertia and the equilibrium magnetization of a ferrofluid. Journal of Magnetism and Magnetic Materials, v. 434, p. 91-99, 2017.

- Gontijo, R. G.; CUNHA, F. R. . Numerical simulations of magnetic suspensions with hydrodynamic and dipole-dipole magnetic interactions. PHYSICS OF FLUIDS, v. 29, p. 062004, 2017.

- GONTIJO, R.G.; Malvar, S. . Microstructural transition in an ordered set of magnetic spheres immersed in a carrier liquid. Mechanics Research Communications, v. 83, p. 12-17, 2017.

- GONTIJO, R.G.. Heat transfer increase for a laminar pipe flow of a magnetic fluid subjected to constant heat flux: an initial theoretical approach. MECHANICS RESEARCH COMMUNICATIONS, v. 91, p. 27-32, 2018.

- DE CARVALHO, D D ; GONTIJO, R G . Reconstructing a continuous magnetization field based on local vorticity cells, CFD and Langevin dynamics: a new numerical scheme. JOURNAL OF MAGNETISM AND MAGNETIC MATERIALS, v. 514, p. 167135, 2020.

- GUIMARÃES, A. B. ; CUNHA, F. R. ; Gontijo, R. G. . The influence of hydrodynamic effects on the complex susceptibility response of magnetic fluids undergoing oscillatory fields: New insights for magnetic hyperthermia. PHYSICS OF FLUIDS, v. 32, p. 012008-012008-17, 2020.

Prof. Rafael Gabler Gontijo, PhD
