#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#CRTBP propagator, with thrust given

function CRTBP_prop_EP_deriv(t, state, MU, DU, TU, Isp, control, time_direction)
    #pull out from state vector:
    x = state[1];
    y = state[2];
    z = state[3];
    xdot = state[4];
    ydot = state[5];
    zdot = state[6];

    if length(state) == 7
        m = state[7];
    else
        m = 1000.0; #default mass (kg)
    end

    #calculate some stuff:
    r1 = sqrt((x+MU)^2 + y^2 + z^2);
    r2 = sqrt((x+MU-1)^2 + y^2 + z^2);

    #calculate r1^3 and r2^3 once each
    r1_3 = r1.^3;
    r2_3 = r2.^3;

    #thrust acceleration magnitude (DU/TU^2):
    T_mag = norm(control)/m/1e3 * TU^2 / DU; # (kg*m/s^2) -> (km/s^2) -> (DU/TU^2) acceleration

    #thrust acceleration vector:
    if norm(control) == 0
        T = control
    else
        T = control ./ norm(control) * T_mag
    end

    g0 = 9.81; #m/s
    mdot = - time_direction * norm(control)/(Isp*g0) * TU; #mass flow rate (kg/TU)

    #allows propagation either forward or backward in the CRTBP:
    omega = time_direction; #+1 for forwards, -1 for backwards

    #calculate accelerations:
    xddot = -(1-MU)*(x+MU)/r1_3 - MU*(x-1+MU)/r2_3 + 2*omega*ydot + x + T[1];
    yddot = -(1-MU)*y/r1_3 - MU*y/r2_3 - 2*omega*xdot + y + T[2];
    zddot = -(1-MU)*z/r1_3 - MU*z/r2_3 + T[3];

    #put back into dstate vector:
    if length(state) == 7
        dstate = [xdot; ydot; zdot; xddot; yddot; zddot; mdot];
    else
        dstate = [xdot; ydot; zdot; xddot; yddot; zddot]; #constant mass
    end

    #Outputs:
    dstate
end



#In-place update of derivatives:
function CRTBP_prop_EP_deriv!(t, state, params, dstate)

    (MU, DU, TU, Isp, control, time_direction) = params

    #pull out from state vector:
    x = state[1];
    y = state[2];
    z = state[3];
    xdot = state[4];
    ydot = state[5];
    zdot = state[6];

    if length(state) == 7
        m = state[7];
    else
        m = 1000.0; #default mass (kg)
    end

    #calculate some stuff:
    r1 = sqrt((x+MU)^2 + y^2 + z^2);
    r2 = sqrt((x+MU-1)^2 + y^2 + z^2);

    #calculate r1^3 and r2^3 once each
    r1_3 = r1.^3;
    r2_3 = r2.^3;

    #thrust acceleration magnitude (DU/TU^2):
    T_mag = norm(control)/m/1e3 * TU^2 / DU; # (kg*m/s^2) -> (km/s^2) -> (DU/TU^2) acceleration

    #thrust acceleration vector:
    if norm(control) == 0
        T = control
    else
        T = control ./ norm(control) * T_mag
    end

    g0 = 9.81; #m/s
    mdot = - time_direction * norm(control)/(Isp*g0) * TU; #mass flow rate (kg/TU)

    #allows propagation either forward or backward in the CRTBP:
    omega = time_direction; #+1 for forwards, -1 for backwards

    #calculate accelerations:
    dstate[1] = xdot
    dstate[2] = ydot
    dstate[3] = zdot
    dstate[4] = -(1-MU)*(x+MU)/r1_3 - MU*(x-1+MU)/r2_3 + 2*omega*ydot + x + T[1]
    dstate[5] = -(1-MU)*y/r1_3 - MU*y/r2_3 - 2*omega*xdot + y + T[2];
    dstate[6] = -(1-MU)*z/r1_3 - MU*z/r2_3 + T[3];

    #put back into dstate vector:
    if length(state) == 7
        dstate[7] = mdot
    end

    #Outputs:
    nothing
end


#for new version of DifferentialEquations, with control computed as a function
#of interpolated costates:
function CRTBP_prop_EP_NNControl_deriv!(dstate, state, params, t)

    (MU, DU, TU, Isp, time_direction, λv_itp, thrustLimit, p, rho) = params


    #λv_itp is an interpolation object for the control. It is prepared like this:
    # using Interpolations
    # λv_itp = scale(interpolate(λv_all', (BSpline(Cubic(Natural())), NoInterp()), OnGrid()), t_TU, 1:3)
    # (from https://github.com/JuliaMath/Interpolations.jl)
    #where t_sec has to be an evenly-spaced time grid.

    #interpolate control from given history:
    λv = [λv_itp[t, 1], λv_itp[t, 2], λv_itp[t, 3]]

    accelLimit = thrustLimit / mass / 1e3 * TU^2 / DU #(kg*m/s^2) -> (km/s^2) -> (DU/TU^2) acceleration

    #Thrust magnitude from control law:
    if p == 0
      #Force thrust to always be on, in the optimal direction. This limits how big
      #of a search we need to make for the neural net stuff.
      umag = accelLimit #N

    elseif p == 1 #minimum mass problem
        g = norm(λv) - 1
        umag = 1/2 * (1 + tanh( g / (2*rho))) * accelLimit

    elseif p > 1 #p is in the range (1, 2]
        umag = (1/p * norm(λv))^(1 / (p-1));

        if umag > accelLimit
          umag = accelLimit
        end
    else
      error("Invalid value of p!")
    end

    #Control vector, scaled properly:
    # control_accel = -umag * λv ./ norm(λv) / mass / 1e3 * TU^2 / DU #(kg*m/s^2) -> (km/s^2) -> (DU/TU^2) acceleration
    control_accel = -umag * λv ./ norm(λv) #DU/TU^2

    if isnan(control_accel[1])
        control_accel[:] = 0.
    end

    #pull out from state vector:
    x = state[1];
    y = state[2];
    z = state[3];
    xdot = state[4];
    ydot = state[5];
    zdot = state[6];

    if length(state) == 7
        m = state[7];
    else
        m = 1000.0; #default mass (kg)
    end

    #calculate some stuff:
    r1 = sqrt((x+MU)^2 + y^2 + z^2);
    r2 = sqrt((x+MU-1)^2 + y^2 + z^2);

    #calculate r1^3 and r2^3 once each
    r1_3 = r1.^3;
    r2_3 = r2.^3;

    g0 = 9.81; #m/s
    mdot = - time_direction * norm(umag)/(Isp*g0) * TU; #mass flow rate (kg/TU)

    #allows propagation either forward or backward in the CRTBP:
    omega = time_direction; #+1 for forwards, -1 for backwards

    #calculate accelerations:
    dstate[1] = xdot
    dstate[2] = ydot
    dstate[3] = zdot
    dstate[4] = -(1-MU)*(x+MU)/r1_3 - MU*(x-1+MU)/r2_3 + 2*omega*ydot + x + control_accel[1]
    dstate[5] = -(1-MU)*y/r1_3 - MU*y/r2_3 - 2*omega*xdot + y + control_accel[2];
    dstate[6] = -(1-MU)*z/r1_3 - MU*z/r2_3 + control_accel[3];

    #if mass is variable:
    if length(state) == 7
        dstate[7] = mdot
    end

    #Outputs:
    nothing
end
