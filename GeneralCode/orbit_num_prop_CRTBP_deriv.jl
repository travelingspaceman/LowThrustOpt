#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Differential equations for CRTBP. Propagate ballistically
#
#Rotating frame


function orbit_num_prop_CRTBP_deriv!(dstate, state, params, t)
    #For use with "DifferentialEquations.jl"

    #pull out parameters:
    MU = params[1]

    #allows propagation either forward or backward in the CRTBP:
    time_direction = params[2] #+1 for forwards, -1 for backwards

    # States
    x   = state[ 1 ];
    y   = state[ 2 ];
    z   = state[ 3 ];
    Vx  = state[ 4 ];
    Vy  = state[ 5 ];
    Vz  = state[ 6 ];

    #ballistic:
    Tx = 0;
    Ty = 0;
    Tz = 0;

    # Compute r1 = Earth-spacecraft
    r1  = ((x+MU).*(x+MU) + y.*y + z.*z).^(1/2);
    r1_3 = r1.^3;

    # Compute r2 = Moon-spacecraft
    r2  = ((x-(1-MU)).*(x-(1-MU)) + y.*y + z.*z).^(1/2);
    r2_3 = r2.^3;

    # Update state derivatives in-place:
    dstate[1] = Vx
    dstate[2] = Vy
    dstate[3] = Vz
    dstate[4] = -(1-MU)*(x+MU)/r1_3 - MU*(x-1+MU)/r2_3 + 2*time_direction*Vy + x + Tx
    dstate[5] = -(1-MU)*y/r1_3 - MU*y/r2_3 - 2*time_direction*Vx + y + Ty
    dstate[6] = -(1-MU)*z/r1_3 - MU*z/r2_3 + Tz
end
