# Dynamical Topography

##Input:

To execute the simulation you must open the Julia's REPL and load the main code:

	julia > include("main.jl")

Once the main code are loaded just input:

	julia > Stokes(n,nSteps)

where ***n*** (blank: 10) is the dimension of the mesh (n)x(n) with Q2Q1-elements. Each element contains ***ppe*** particles for each one and ***nSteps*** (blank: 20) are the time steps.

Alternatively:

	julia > Stokes(n,nSteps,air)

where air adding a layer of air on the top of the domain.

##Output

The output of the two case are one data files on `./results/file.dat` and one figure of the last map of the particles on `./results/figs/fig.png`.

## Recover the Data

It is possible to recover the final height-plots using:

	julia > recoverdata(n,case,nSteps)

where ***n*** is the dimension of the mesh (n)x(n) (specified in the name of the file) and the ***case*** is 1 for the *Fixed* case and *Air Layer* case.

## Modifying the Code
You can change the following data modifying it in Stokes.jl source code and reloaded the main code *main.jl*:

Warning with these variations, the code is not plenty tested.

- The domain ***x1,x2,y1,y2*** (default: 40x40, mandatory square)
- The number of particles per element ***ppe*** (default: 30)
- The position and the radius of the sphere ***pto, radius*** (default: (20,10), 1.75)
- The densities ***rho*** and viscosities ***visc*** (default: (1420,1150), (50,69000))
- The radius of tolerance ***sradius*** must changes depends on ***n***.
