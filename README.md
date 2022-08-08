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

__1. What does it simulate?__

Sims was designed to perform simultaneous numerical experiments of N magnetic particles immersed in a viscous liquid during sedimentation or shear induced motion. The particles can be either neutrally buoyant or denser than the base fluid. They can undergo a steady state or either an unsteady oscillatory shear field. They can also be in the presence or in the absence of an applied magnetic field. When in the presence of an applied field, this field can be either steady-state or time dependent. In the case of time dependent field excitation the user can choose an oscillatory 1D field, a rotating 2D field, a non-linear 1D Duffing excitation or even a 1D oscillatory beating pattern excitation. The user can also sets other options, such as: compute or neglect Brownian effects, change the initial configuration of the particles to a set of ordered spheres or to a random dispersion of spheres inside a spherical aggregate-like structure, mix magnetic with non-magnetic particles and compute only velocity fluctuations. 

__2. The configuration file (entrada.dat)__

The user defines what scenario should be simulated by changing the data displayed in the configuration file 'entrada.dat'. In the configuration file the first thing is to define the values of important logical variables (true or false) by basically answering some questions regarding what is going to be simulated. These answers will be used in conditional structures all along the code. In the configuration file the user can also sets numerical data regarding the number of particles and the number of independent realizations, the volume fraction of particles, the aspect lenght of the simulation box and the number of physical and reciprocal lattice cells used to compute long range interactions, which are relevant only when we compute hydrodynamic interactions and dipolar magnetic torques in concentrated regimes where the volume fraction of particles is higher than 10%. The number 125 for this option should work fine. The user can also set the frequencies of the applied shear rate and magnetic field, when it is the case. These frequencies are treated in a nondimensional form, using characteristic scales defined based on the logical answers at the beginning of the configuration file. 

Informations regarding the intensity of the applied field, magnitude of magnetic moments, size of the particles, viscosity of the surrounding liquid and local temperature are expressed in terms of non-dimensional physical parameters, namely:

St = Stokes number

Str = Rotational Stokes number

Pe = Péclet number

alpha = dimensionless magnetic field

lambda = dimensionless intensity of the magnitude of the magnetic moments of the particles

Detailed information regarding these physical parameters can be found on the References at the end of this readme file.

__3. Data file output__

In the configuration file the user may choose to record the positions, velocities and/or orientation of the dipole moments of the particles in external data files. If the option "RECORD POSITION IN FILE" is set as TRUE the solver will create several files named "posicao  x.plt" where "x" is an index that represents the number of that particular numerical experiment. These files contain the position in each direction (X,Y,Z) of each particle (each line) for several time-steps (separated by zones). If the option "RECORD DIPOLE IN FILE" is enabled the orientations of each magnetic moment will be included in the "posicao  x.plt" files. Only these files are necessary to make animations of the motion of the particles using post-processing softwares such as Tecplot and Ovito.

If the option "RECORD VELOCITY IN FILE" is set to be true, the solver will create several files named "velocidade  x.plt" where "x" is an index that represents the number of that particular numerical experiment. These files are not necessary for particle animation, but are required to perform post-processing statistical analysis. So, if the option "STATISTICAL ANALYSIS" is enabled, than "RECORD POSITION IN FILE", "RECORD VELOCITY IN FILE" and "RECORD DIPOLE IN FILE" must all be true.

The option "CALCULATE MAGNETIZATION" is also very important here. If this option is set to be TRUE, then the solver will automatically create a data file named mag_tempo.plt, which will be updated during the simulation. This file contains the average orientation (magnetization) of all the particles in all simultaneous numerical experiments in each direction (Mx, My, Mz), the field excitation components and its derivatives as a function of time.

__4. Magnetic hyperthermia specificities__

In the problem of Magnetic Hyperthermia we are interested in computing the response of a collection of magnetic particles to the application of a time-dependent field. The link between excitation (H) and response (M) is directly related to the capacity of the multibody system in dissipating the field's energy in the form of heat. Therefore, the M x H relation is an important relation, that can be directly extracted from the "mag_tempo.plt" file. Moreover, the integral of M.dH, which represents the area in a hysteresis loop is also a very important information. This area is presented in the terminal after the end of the simulation in the form of a short MHT REPORT. We have also investigated these multibody systems using dynamical system tools. For this purpose we have implemented the possibility of an automatic frequency sweep in order to build some bifurcations diagrams of the system.


__5. References__

For more information regarding the mathematical formulation used in SIMS, examples and a better physical interpretation of the dimensionless parameters mentioned above, please consult the following references.

[1] GONTIJO, R.G.; CUNHA, F.R. . Dynamic numerical simulations of magnetically interacting suspensions in creeping flow. Powder Technology (Print), v. 279, p. 146-165, 2015.

[2] GONTIJO, R.G.; Malvar, S. ; CUNHA, F.R. . Magnetic particulate suspensions from the perspective of a dynamical system. Powder Technology (Print), v. 297, p. 165-182, 2016.

[3] GONTIJO, R.G.. A numerical perspective on the relation between particle rotational inertia and the equilibrium magnetization of a ferrofluid. Journal of Magnetism and Magnetic Materials, v. 434, p. 91-99, 2017.

[4] Gontijo, R. G.; CUNHA, F. R. . Numerical simulations of magnetic suspensions with hydrodynamic and dipole-dipole magnetic interactions. PHYSICS OF FLUIDS, v. 29, p. 062004, 2017.

[5] GONTIJO, R.G.; Malvar, S. . Microstructural transition in an ordered set of magnetic spheres immersed in a carrier liquid. Mechanics Research Communications, v. 83, p. 12-17, 2017.

[6] GONTIJO, R.G.. Heat transfer increase for a laminar pipe flow of a magnetic fluid subjected to constant heat flux: an initial theoretical approach. MECHANICS RESEARCH COMMUNICATIONS, v. 91, p. 27-32, 2018.

[7] DE CARVALHO, D D ; GONTIJO, R G . Reconstructing a continuous magnetization field based on local vorticity cells, CFD and Langevin dynamics: a new numerical scheme. JOURNAL OF MAGNETISM AND MAGNETIC MATERIALS, v. 514, p. 167135, 2020.

[8] GUIMARÃES, A. B. ; CUNHA, F. R. ; Gontijo, R. G. . The influence of hydrodynamic effects on the complex susceptibility response of magnetic fluids undergoing oscillatory fields: New insights for magnetic hyperthermia. PHYSICS OF FLUIDS, v. 32, p. 012008-012008-17, 2020.

Prof. Rafael Gabler Gontijo, PhD
