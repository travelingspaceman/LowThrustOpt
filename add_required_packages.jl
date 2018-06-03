#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Commonly used 3rd party packages
#Install all in one swoop.

#Update all packages -- NOTE this is dangerous, as Julia is still a young
#language. Updates frequently break things:
# Pkg.update()

Pkg.add("DifferentialEquations")
# Pkg.add("Gurobi") #Optional, for direct method, requires license
Pkg.add("Plots")
Pkg.add("PlotlyJS")
Pkg.add("Ipopt") #Open-source alternative to Gurobi
# Pkg.add("JPLEphemeris") #for interplanetary problems, not CRTBP
# Pkg.add("MATLAB") #Optional, used for plotting, requires license
Pkg.add("JuMP")

#Trigger pre-compilation:
using DifferentialEquations
# using Gurobi
using Plots
using PlotlyJS
using Ipopt
# using JPLEphemeris
# using MATLAB
using JuMP
