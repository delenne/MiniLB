# MiniLB
Minimalist LBM Code

This code compute a Poiseuille flow into a pipe using Lattice Boltzmann Method (LBM) with D2Q9 scheme and BJK operator.

## Parameters
All physical parameters are included at the begining of the main file.

### Dimentions of the pipe 
length: lx 
size (height): ly 

### Pressure gradient at left and right boder
rhoin: incomming densities
rhoout: outgoing densities

For the pipe wall, bounce back conditions (no-slip) are used for top and bottom boundaries.

### Important note
For the sake of clarity only basic boundary conditions have been integreted.
To avoid instabilities for higher pressure drop, specific boundary conditions
should be set for the 4 corners of the domain.

### Relaxation time
tau

### Output results
period: save result files every 'period' 

## Compilation

In a terminal
g++ -O3 ./main.cpp -o MiniLB

to run 
./MiniLB

## Visualisation

The code outputs periodically two type of files:

### Pressure and Velocity fields
output_P_U_Fields_???.vtk

These file can be open using Paraview (https://www.paraview.org)

### Velocity profile
output_U_Profile_???.txt
This file contains the velocity profile in section perpendicular to the fluid flow located at in the midle of the pipe. 


