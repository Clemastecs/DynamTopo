

"""
	This script load the workspace of Stokes' model for dynamical topography.

"""

#workspace() # Clean workspace memory

## Packages ##
##############
# Use Pkg.add("name_pkg") to install if not
# Use Pkg.update() to update the packages
# We need explicit call for the import pkgs. e.i. PyPlot.function()

using Grid # Grid Interpolation, instructions: https://github.com/timholy/Grid.jl
import PyPlot # Plots, install python-mathplotlib first


## Libraries ##
###############
# (Calling our scripts located in the same directory of main.jl)

include("BC_freeSlip.jl")
include("CreateVelNod.jl")
include("CreatePresNod.jl")
include("CreMat.jl")
include("EleMat.jl")
include("Source.jl")
include("plotLayer.jl")
include("shapeFunc.jl")
include("plotPar.jl")
include("setparticles.jl")
include("Quadrature.jl")
include("updateParticles!.jl")
include("remesh!.jl")
include("compStress.jl")
include("plotFix.jl")

include("Stokes.jl")

include("recoverdata.jl")

## Body ##
##########

air=bool(true)


print_with_color(:blue,"====>\n")
print_with_color(:blue," (i) Input Stokes(n,nSteps) to run a nSteps time steps simulation in a domain of [0,40]x[0,40] with (n)x(n) Q2Q1-elements with aprox. 20 particles for each one\n")
print("\n")
print_with_color(:blue," (ii) Input Stokes(n,nStep,air) to run a simulation in a domain of [0,40]x[0,40] with (n)x(n+1*n/4) Q2Q1-elements with aprox. 20 particles for each one, adding a layer of air\n")
print("\n")
print_with_color(:blue," OUTPUT: One data files on [./results/file.dat] and one figure of the last map of the particles on [./results/figs/fig.png] . \n")
print("\n")
print_with_color(:blue," n default = 10\n")
print_with_color(:blue," nSteps default = 20\n")
print_with_color(:blue," air default = false\n")
print_with_color(:blue,"<====\n")
print("\n")
print_with_color(:blue,"(See the README.md file for more details.)\n")
