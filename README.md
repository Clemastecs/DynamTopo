# Dynamical Topography

A workspace to compute a Stokes' model for dynamical topography, combining FEM/Marker-in-Cell.

##Input:

To execute the simulation you must open the Julia's REPL and load the main code:

	julia > include("main.jl")

Once the main code is loaded just input:

	julia > Stokes(n,nSteps)

where ***n*** is the dimension of the mesh (n)x(n) with Q2Q1-elements. Each element contains ***ppe*** particles and ***nSteps*** are the time steps.

Alternatively:

	julia > Stokes(n,nSteps,air)

where ***air*** adds a layer of "sticky" air on the top of the domain.

##Output

The outputs of the two cases are one data files on `./results/file.dat` and one figure of the last map of the particles on `./results/figs/fig.png`.

## Recover the Data

It is possible to recover the final height-plots using:

	julia > recoverdata(n,case,nSteps)

where ***n*** is the dimension of the mesh (n)x(n) (specified in the name of the file) and ***case*** is 1 for the *Fixed* case and 2 for the *Air Layer* case.

## Annotations

This simulation needs the large number of elements to works fine. It is mandatory uses a square domain.

## Future releases

	- Re-meshing the domain with a non-regular mesh with more definition in a neighbourhood of the boundary between the "sticky" air layer and the mantle (It needs a non-regular grid's interpolation).

## Modifying the Code

You can change the following data modifying it in Stokes.jl source code and reloading the main code *main.jl*:

Warning with these variations, the code is not plenty tested.

- The number of particles per element ***ppe*** (default: 30)
- The position and the radius of the sphere ***pto, radius***
- The densities ***rho*** and the viscosities ***visc*** (default: (1420,1150), (50,69000))
- The radius of tolerance ***sradius*** must change depending on ***n***.
