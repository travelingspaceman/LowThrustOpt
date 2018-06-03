#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#
#This module, "LowThrustOpt", is designed for low-thrust spacecraft trajectory
#optimization, especially in the Earth-Moon circular restricted 3-body problem
#dynamics. The algorithms implemented here are also valid in any other
#dynamical environment.
#
#Direct multiple shooting
#Indirect multiple shooting



#To add the LowThrustOpt directory to the path for a single session:
# push!(LOAD_PATH, "C:/.../path/to/containing/folder")


module LowThrustOpt

#Earth-Moon system:
const global MU = 0.012150585609624037 # MU = mu_moon/(mu_moon + mu_planet)
const global DU = 384747.96285603708 #Earth-Moon distance (km) (km per distance unit)
const global TU = 375699.81732246041 #time units (seconds per time unit)
const global r_moon = 1737. #km
const global r_earth = 6378. #km
const global day = 86400. #seconds
##Calculate mu_moon based on MU for consistency with other code:
const global mu_planet = 398600.4415
const global mu_moon = (MU * mu_planet) / (1 - MU)

#External packages required:
using DifferentialEquations
using Ipopt
using JuMP
using Interpolations
using ForwardDiff
using GeneralCode
using Plots

export multiShoot_CRTBP_direct
export plotTrajPlotly_direct
export meshRefine_direct
export multiShoot_CRTBP_indirect
export plotTrajPlotly_indirect
export controlLaw_cart
export CRTBP_stateCostate_deriv!
export CRTBP_prop_EP_NNControl_deriv!
export jacobiConstant
export interpInitialStates
export find_τ
export densify
export reduceFuel_indirect
export addTimeFinal


include("multiShoot_CRTBP_direct.jl")
include("multiShoot_CRTBP_indirect.jl")
include("CRTBP_stateCostate_deriv.jl")
include("CRTBP_prop_EP_deriv.jl")
include("HelperFunctions.jl")


end
