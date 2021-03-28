# MiniLB
Minimalist LBM Code

This code compute a Poiseuille flow into a pipe using Lattice Boltzmann Method (LBM) with D2Q9 scheme and BJK operator.
All physical parameters are included at the begining of the main file.

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

