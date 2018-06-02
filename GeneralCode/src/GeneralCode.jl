#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Module for a collection of generally useful functions.

#To add the GeneralCode directory to the path for a single session:
# push!(LOAD_PATH, "C:/.../path/to/containing/folder")

#Declare module name:
module GeneralCode

    #List of functions exported by this module:
    export meeus_ephem
    export norm_many
    export orbit_num_prop_deriv
    export orbit_num_prop_deriv!
    export rv2coes
    export LinInterp
    export ode4
    export ode5
    export ode7
    export ode78
    export ode78_events
    export ode7_8
    export ode7_8!
    export orbit_num_prop_CRTBP_deriv
    export orbit_num_prop_CRTBP_deriv!
    export rv2mee
    export rv2mee_r
    export mee2rv
    export vector_rotate
    export skewSymmetric
    export synodic2inertial
    export inertial2synodic
    export sphere
    export TwoBody_prop_EP_deriv
    export TwoBody_prop_EP_deriv!
    export TwoBody_prop_EP_precomputedControl_deriv!
    export TwoBody_prop_EP_NNControl_deriv!
    export coes2rv3
    export twoBody_stateCostate_deriv
    export twoBody_stateCostate_deriv!
    export twoBody_stateCostate_deriv_STM!
    export twoBody_stateCostate_mass_deriv
    export twoBody_stateCostate_mass_deriv!
    export long_short_way
    export lambert
    export covarianceEmpirical
    export cart2RaDec
    export RaDec2cart
    export orbit_analy_prop


    #List of Julia scripts included in this module:
    include("meeus.jl")
    include("norm_many.jl")
    include("orbit_num_prop_deriv.jl")
    include("rv2coes.jl")
    include("LinInterp.jl")
    include("ode.jl")
    include("orbit_num_prop_CRTBP_deriv.jl")
    include("rv2mee.jl")
    include("mee2rv.jl")
    include("vector_rotate.jl")
    include("skewSymmetric.jl")
    include("CRTBP_synodicInertialConvert.jl")
    include("sphere.jl")
    include("TwoBody_prop_EP_deriv.jl")
    include("coes2rv3.jl")
    include("twoBody_stateCostate_deriv.jl")
    include("twoBody_stateCostate_mass_deriv.jl")
    include("long_short_way.jl")
    include("lambert.jl")
    include("covarianceEmpirical.jl")
    include("cart2RADEC.jl")
    include("orbit_analy_prop.jl")


end
