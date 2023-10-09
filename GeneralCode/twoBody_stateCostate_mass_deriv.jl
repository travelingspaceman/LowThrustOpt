#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Derivatives for states & costates -- for interplanetary low-thrust shooting
#


#for NEW version of 'DifferentialEquations.jl'
function twoBody_stateCostate_mass_deriv!(dstate_costate, state_costate, params, t)
    #for use with "DifferentialEquations.jl" or "ode7_8!"

    #pull out parameters
    (mu_C, Isp, p, PowerSystem, thrustFactor, ρ) = params

    V = state_costate[4:6]
    mass = state_costate[7]

    lambda_v = state_costate[11:13];

    #power available at time t for thrusting (more accurate model)
    AU = 149597870.700 #km
    thrustLimit_reduced = PowerSystem( t / (86400 * 365.25), norm(state_costate[1:3] / AU), thrustFactor )

    accelLimit = thrustLimit_reduced / mass / 1e3 # km/s2

    #Thrust magnitude from control law:
    if p == 0
        #Force thrust to always be on, in the optimal direction. This limits how big
        #of a search we need to make for the neural net stuff.
        umag = accelLimit #km/s^2

    elseif p == 1 #minimum mass problem
        #Smooth the discontinuity with hyperbolic tangent:
        g = norm(lambda_v) - 1
        umag = 1/2 * (accelLimit + accelLimit*tanh(g/ρ) )

    elseif p > 1 #p is in the range (1, 2]
        umag = (1/p * norm(lambda_v))^(1 / (p-1));

        if umag > accelLimit
            umag = accelLimit
        end
    else
        error("Invalid value of p!")
    end

    #Control vector, scaled properly:
    control_accel = -umag * lambda_v ./ norm(lambda_v) #km/s^2
    if isnan(control_accel[1])
        control_accel = zeros(3)
    end

    #mass flow rate:
    g0 = 9.80665 #m/s2
    mdot = -(umag * mass * 1e3) / (Isp * g0) #kg/s

    dstate_costate[1:3] = state_costate[4:6]
    dstate_costate[4:6] = -mu_C/norm(state_costate[1:3])^3 * state_costate[1:3] + control_accel
    dstate_costate[7] = mdot

    L4 = state_costate[11]
    L5 = state_costate[12]
    L6 = state_costate[13]

    X1 = state_costate[1]
    X2 = state_costate[2]
    X3 = state_costate[3]

    rmag5 = norm(state_costate[1:3])^5 #compute once and use 3 times
    dstate_costate[8] = -(mu_C*(2*L4*X1^2 + 3*L5*X1*X2 + 3*L6*X1*X3 - L4*X2^2 - L4*X3^2))/rmag5
    dstate_costate[9] = -(mu_C*(- L5*X1^2 + 3*L4*X1*X2 + 2*L5*X2^2 + 3*L6*X2*X3 - L5*X3^2))/rmag5
    dstate_costate[10] = -(mu_C*(- L6*X1^2 + 3*L4*X1*X3 - L6*X2^2 + 3*L5*X2*X3 + 2*L6*X3^2))/rmag5
    dstate_costate[11:13] = -state_costate[8:10]
    dstate_costate[14] = dot(state_costate[11:13], control_accel) / (1e3*mass)

end
