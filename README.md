# Lattice Boltzmann Method

This code simulates D2Q9 (two dimension with nine velocities) Lattice Boltzmann fluid flow on two dimensional square lattice in the presence of circular/rectangular obstacle. Initial drift is set along positive x direction. Colorplot of vorticity is shown for both circular and rectangular obstacles. Karman vortices are 
being observed in the attached plots. We have used reflecting/bounce-back boundary conditions while encountering with solid surfaces. In the attached plots
periodic boundary conditions are applied along both x and y directions. The code also has the facility of switching to reflecting-y boundary conditions. Using
gnuplot pipe live visualisation of different system variables (velocity and vorticlty) are possible.

![vorticity_circular](https://github.com/dipanjan308/lattice_boltzmann_method_d2q9/assets/61291745/a35104df-2a12-4406-a452-0c6003663489)
![vorticity_rectangular](https://github.com/dipanjan308/lattice_boltzmann_method_d2q9/assets/61291745/dbf6821b-8aa3-4114-b8b8-50143358a505)

