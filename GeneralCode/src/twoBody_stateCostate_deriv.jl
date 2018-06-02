#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Derivatives for states & costates -- for Earth-Mars low-thrust
#single-shooting
#
#Constant mass, constant thrust


#For ode78()
function twoBody_stateCostate_deriv(t, state_costate, mu_C, mass, p, PowerSystem)

    V = state_costate[4:6];

    lambda_r = state_costate[7:9];
    lambda_v = state_costate[10:12];


    #power available at time t for thrusting (more accurate model)
    AU = 149597870.700 #km
    thrustLimit_reduced = PowerSystem( t / (86400 * 365.25), norm(state_costate[1:3] / AU), thrustFactor )

    #Thrust magnitude from control law:
    if p == 0
        #Force thrust to always be on, in the optimal direction. This limits how big
        #of a search we need to make for the neural net stuff.
        umag = thrustLimit_reduced #N

    elseif p == 1 #minimum mass problem
        if norm(lambda_v) >= 1
            umag = thrustLimit_reduced
        else
            umag = 0.
        end

        #Smooth the discontinuity with hyperbolic tangent:
        g = norm(lambda_v) - 1
        umag = 1/2 * (1 + tanh(g/(2*rho) ) ) * thrustLimit_reduced


    elseif p > 1 #p is in the range (1, 2]
        umag = (1/p * norm(lambda_v))^(1 / (p-1));

        if umag > thrustLimit_reduced
            umag = thrustLimit_reduced
        end
    else
        error("Invalid value of p!")
    end

    #Control vector, scaled properly:
    control_accel = -umag * lambda_v ./ norm(lambda_v) / mass / 1e3 #(kg*m/s^2) -> (km/s^2) acceleration
    if isnan(control_accel[1])
        control_accel = zeros(3)
    end

    r1 = state_costate[1];
    r2 = state_costate[2];
    r3 = state_costate[3];
    rmag = norm(state_costate[1:3])

    L4 = lambda_v[1];
    L5 = lambda_v[2];
    L6 = lambda_v[3];

    #Output:
    dstate_costate = zeros(12)
    dstate_costate[1:3] = state_costate[4:6]
    dstate_costate[4:6] = -mu_C/norm(state_costate[1:3])^3 * state_costate[1:3] + control_accel

    dstate_costate[7]    = (mu_C.*(2.*L4.*r1.^2 + 3.*L5.*r1.*r2 + 3.*L6.*r1.*r3 - L4.*r2.^2 - L4.*r3.^2))./(rmag).^5
    dstate_costate[8]    = (mu_C.*(- L5.*r1.^2 + 3.*L4.*r1.*r2 + 2.*L5.*r2.^2 + 3.*L6.*r2.*r3 - L5.*r3.^2))./(rmag).^5
    dstate_costate[9]    = (mu_C.*(- L6.*r1.^2 + 3.*L4.*r1.*r3 - L6.*r2.^2 + 3.*L5.*r2.*r3 + 2.*L6.*r3.^2))./(rmag).^5
    dstate_costate[10:12] = state_costate[7:9]

    dstate_costate
end



#For NEW version of 'DifferentialEquations.jl'
function twoBody_stateCostate_deriv!(dstate_costate, state_costate, params, t)

    #pull out parameters
    (mu_C, mass, p, PowerSystem, thrustFactor, rho) = params

    V = state_costate[4:6];

    lambda_r = state_costate[7:9];
    lambda_v = state_costate[10:12];

    #power available at time t for thrusting (more accurate model)
    AU = 149597870.700 #km
    thrustLimit_reduced = PowerSystem( t / (86400 * 365.25), norm(state_costate[1:3] / AU), thrustFactor )

    #Thrust magnitude from control law:
    if p == 0
        #Force thrust to always be on, in the optimal direction. This limits how big
        #of a search we need to make for the neural net stuff.
        umag = thrustLimit_reduced #N

    elseif p == 1 #minimum mass problem

        #Smooth the discontinuity with hyperbolic tangent:
        g = norm(lambda_v) - 1
        umag = 1/2 * (1 + tanh(g/ (2*rho) ) ) * thrustLimit_reduced

    elseif p > 1 #p is in the range (1, 2]
        umag = (1/p * norm(lambda_v))^(1 / (p-1));

        if umag > thrustLimit_reduced
            umag = thrustLimit_reduced
        end
    else
        error("Invalid value of p!")
    end



    #Control vector, scaled properly:
    control_accel = -umag * lambda_v ./ norm(lambda_v) / mass / 1e3 #(kg*m/s^2) -> (km/s^2) acceleration
    if isnan(control_accel[1])
        control_accel = zeros(3)
    end

    r1 = state_costate[1];
    r2 = state_costate[2];
    r3 = state_costate[3];
    rmag = norm(state_costate[1:3])

    L4 = lambda_v[1];
    L5 = lambda_v[2];
    L6 = lambda_v[3];

    #Output in place:
    dstate_costate[1:3] = state_costate[4:6]
    dstate_costate[4:6] = -mu_C/norm(state_costate[1:3])^3 * state_costate[1:3] + control_accel

    dstate_costate[7] = (mu_C.*(2.*L4.*r1.^2 + 3.*L5.*r1.*r2 + 3.*L6.*r1.*r3 - L4.*r2.^2 - L4.*r3.^2))./(rmag).^5
    dstate_costate[8] = (mu_C.*(- L5.*r1.^2 + 3.*L4.*r1.*r2 + 2.*L5.*r2.^2 + 3.*L6.*r2.*r3 - L5.*r3.^2))./(rmag).^5
    dstate_costate[9]    = (mu_C.*(- L6.*r1.^2 + 3.*L4.*r1.*r3 - L6.*r2.^2 + 3.*L5.*r2.*r3 + 2.*L6.*r3.^2))./(rmag).^5
    dstate_costate[10:12] = state_costate[7:9]
end


#for NEW version of 'DifferentialEquations.jl'
function twoBody_stateCostate_deriv_STM!(dstate_costate, state_costate, params, t::Float64)
    #for use with "DifferentialEquations.jl"
    #Uses scaled inputs and outputs!

    #pull out parameters
    (mu_C, mass, p, PowerSystem, thrustFactor, scale_vec, scale_mat_vec) = params
    #scale_mat_vec is equal to: scale_mat[:]


    #un-scale to physical units:
    #have to create a copy variable, otherwise 'state_costate' gets changed in
    #other scopes...
    state_costate_unscaled = copy(state_costate[1:12]) .* scale_vec

    x  = state_costate_unscaled[1]
    y  = state_costate_unscaled[2]
    z  = state_costate_unscaled[3]
    vx = state_costate_unscaled[4]
    vy = state_costate_unscaled[5]
    vz = state_costate_unscaled[6]
    L1 = state_costate_unscaled[7]
    L2 = state_costate_unscaled[8]
    L3 = state_costate_unscaled[9]
    L4 = state_costate_unscaled[10]
    L5 = state_costate_unscaled[11]
    L6 = state_costate_unscaled[12]

    #state transition matrix
    ϕ = reshape(state_costate[13:end] .* scale_mat_vec, 12, 12) #un-scale and reshape

    rmag = norm(state_costate_unscaled[1:3])

    V = state_costate_unscaled[4:6];
    lambda_v = state_costate_unscaled[10:12];
    lambda_v_mag = norm(lambda_v)

    #power available at time t for thrusting (more accurate model)
    AU = 149597870.700 #km
    thrustLimit_reduced = PowerSystem( t / (86400 * 365.25), norm(state_costate_unscaled[1:3] / AU), thrustFactor )

    #Thrust magnitude from control law:
    if p == 1 #minimum fuel mass problem
        if norm(lambda_v) >= 1
            umag = thrustLimit_reduced

            #For the case when thrust is at thrustLimit
            A = vcat(hcat(zeros(3,3), eye(3), zeros(3, 6)),
                [ -(mu_C*(rmag^2 - 3*x*x))/rmag^5, (3*mu_C*x*y)/rmag^5, (3*mu_C*x*z)/rmag^5, 0, 0, 0, 0, 0, 0, -(thrustLimit_reduced*(L5^2 + L6^2))/(1000*mass*lambda_v_mag^3), (L4*L5*thrustLimit_reduced)/(1000*mass*lambda_v_mag^3), (L4*L6*thrustLimit_reduced)/(1000*mass*lambda_v_mag^3)]',
                [ (3*mu_C*y*x)/rmag^5, -(mu_C*(rmag^2 - 3*y*y))/rmag^5, (3*mu_C*y*z)/rmag^5, 0, 0, 0, 0, 0, 0, (L4*L5*thrustLimit_reduced)/(1000*mass*lambda_v_mag^3), -(thrustLimit_reduced*(L4^2 + L6^2))/(1000*mass*lambda_v_mag^3), (L5*L6*thrustLimit_reduced)/(1000*mass*lambda_v_mag^3)]',
                [ (3*mu_C*z*x)/rmag^5, (3*mu_C*z*y)/rmag^5, -(mu_C*(rmag^2 - 3*z*z))/rmag^5, 0, 0, 0, 0, 0, 0, (L4*L6*thrustLimit_reduced)/(1000*mass*lambda_v_mag^3), (L5*L6*thrustLimit_reduced)/(1000*mass*lambda_v_mag^3), -(thrustLimit_reduced*(L4^2 + L5^2))/(1000*mass*lambda_v_mag^3)]',
                [ (3*mu_C*(- 2*L4*x^3 - 4*L5*x^2*y - 4*L6*x^2*z + 3*L4*x*y^2 + 3*L4*x*z^2 + L5*y^3 + L6*y^2*z + L5*y*z^2 + L6*z^3))/rmag^7, (3*mu_C*(L5*x^3 - 4*L4*x^2*y - 4*L5*x*y^2 - 5*L6*x*y*z + L5*x*z^2 + L4*y^3 + L4*y*z^2))/rmag^7, (3*mu_C*(L6*x^3 - 4*L4*x^2*z + L6*x*y^2 - 5*L5*x*y*z - 4*L6*x*z^2 + L4*y^2*z + L4*z^3))/rmag^7, 0, 0, 0, 0, 0, 0, -(mu_C*(- 2*x^2 + y^2 + z^2))/rmag^5, (3*mu_C*x*y)/rmag^5, (3*mu_C*x*z)/rmag^5]',
                [ (3*mu_C*(L5*x^3 - 4*L4*x^2*y - 4*L5*x*y^2 - 5*L6*x*y*z + L5*x*z^2 + L4*y^3 + L4*y*z^2))/rmag^7, (3*mu_C*(L4*x^3 + 3*L5*x^2*y + L6*x^2*z - 4*L4*x*y^2 + L4*x*z^2 - 2*L5*y^3 - 4*L6*y^2*z + 3*L5*y*z^2 + L6*z^3))/rmag^7, (3*mu_C*(L6*x^2*y + L5*x^2*z - 5*L4*x*y*z + L6*y^3 - 4*L5*y^2*z - 4*L6*y*z^2 + L5*z^3))/rmag^7, 0, 0, 0, 0, 0, 0, (3*mu_C*x*y)/rmag^5, -(mu_C*(x^2 - 2*y^2 + z^2))/rmag^5, (3*mu_C*y*z)/rmag^5]',
                [ (3*mu_C*(L6*x^3 - 4*L4*x^2*z + L6*x*y^2 - 5*L5*x*y*z - 4*L6*x*z^2 + L4*y^2*z + L4*z^3))/rmag^7, (3*mu_C*(L6*x^2*y + L5*x^2*z - 5*L4*x*y*z + L6*y^3 - 4*L5*y^2*z - 4*L6*y*z^2 + L5*z^3))/rmag^7, (3*mu_C*(L4*x^3 + L5*x^2*y + 3*L6*x^2*z + L4*x*y^2 - 4*L4*x*z^2 + L5*y^3 + 3*L6*y^2*z - 4*L5*y*z^2 - 2*L6*z^3))/rmag^7, 0, 0, 0, 0, 0, 0, (3*mu_C*x*z)/rmag^5, (3*mu_C*y*z)/rmag^5, -(mu_C*(x^2 + y^2 - 2*z^2))/rmag^5]',
                hcat(zeros(3, 6), eye(3), zeros(3,3)) )
        else
            umag = 0.

            #For the case when there is no control (p=1, switching function off) (need to check that it works, but no reason to believe it wouldn't)
            A = vcat(hcat(zeros(3,3), eye(3), zeros(3, 6)),
                [-(mu_C*(rmag^2 - 3*x*x))/rmag^5, (3*mu_C*x*y)/rmag^5, (3*mu_C*x*z)/rmag^5, 0, 0, 0, 0, 0, 0, 0, 0, 0]',
                [(3*mu_C*y*x)/rmag^5, -(mu_C*(rmag^2 - 3*y*y))/rmag^5, (3*mu_C*y*z)/rmag^5, 0, 0, 0, 0, 0, 0, 0, 0, 0]',
                [(3*mu_C*z*x)/rmag^5, (3*mu_C*z*y)/rmag^5, -(mu_C*(rmag^2 - 3*z*z))/rmag^5, 0, 0, 0, 0, 0, 0, 0, 0, 0]',
                [(3*mu_C*(- 2*L4*x^3 - 4*L5*x^2*y - 4*L6*x^2*z + 3*L4*x*y^2 + 3*L4*x*z^2 + L5*y^3 + L6*y^2*z + L5*y*z^2 + L6*z^3))/rmag^7, (3*mu_C*(L5*x^3 - 4*L4*x^2*y - 4*L5*x*y^2 - 5*L6*x*y*z + L5*x*z^2 + L4*y^3 + L4*y*z^2))/rmag^7, (3*mu_C*(L6*x^3 - 4*L4*x^2*z + L6*x*y^2 - 5*L5*x*y*z - 4*L6*x*z^2 + L4*y^2*z + L4*z^3))/rmag^7, 0, 0, 0, 0, 0, 0, -(mu_C*(- 2*x^2 + y^2 + z^2))/rmag^5, (3*mu_C*x*y)/rmag^5, (3*mu_C*x*z)/rmag^5]',
                [(3*mu_C*(L5*x^3 - 4*L4*x^2*y - 4*L5*x*y^2 - 5*L6*x*y*z + L5*x*z^2 + L4*y^3 + L4*y*z^2))/rmag^7, (3*mu_C*(L4*x^3 + 3*L5*x^2*y + L6*x^2*z - 4*L4*x*y^2 + L4*x*z^2 - 2*L5*y^3 - 4*L6*y^2*z + 3*L5*y*z^2 + L6*z^3))/rmag^7, (3*mu_C*(L6*x^2*y + L5*x^2*z - 5*L4*x*y*z + L6*y^3 - 4*L5*y^2*z - 4*L6*y*z^2 + L5*z^3))/rmag^7, 0, 0, 0, 0, 0, 0, (3*mu_C*x*y)/rmag^5, -(mu_C*(x^2 - 2*y^2 + z^2))/rmag^5, (3*mu_C*y*z)/rmag^5]',
                [(3*mu_C*(L6*x^3 - 4*L4*x^2*z + L6*x*y^2 - 5*L5*x*y*z - 4*L6*x*z^2 + L4*y^2*z + L4*z^3))/rmag^7, (3*mu_C*(L6*x^2*y + L5*x^2*z - 5*L4*x*y*z + L6*y^3 - 4*L5*y^2*z - 4*L6*y*z^2 + L5*z^3))/rmag^7, (3*mu_C*(L4*x^3 + L5*x^2*y + 3*L6*x^2*z + L4*x*y^2 - 4*L4*x*z^2 + L5*y^3 + 3*L6*y^2*z - 4*L5*y*z^2 - 2*L6*z^3))/rmag^7, 0, 0, 0, 0, 0, 0, (3*mu_C*x*z)/rmag^5, (3*mu_C*y*z)/rmag^5, -(mu_C*(x^2 + y^2 - 2*z^2))/rmag^5]',
                hcat(zeros(3, 6), eye(3), zeros(3,3)) )
        end

    elseif p > 1 #p is in the range (1, 2]
        umag = (1/p * norm(lambda_v))^(1 / (p-1));

        #For the case when thrust is between 0 and thrustLimit
        A = vcat(hcat(zeros(3,3), eye(3), zeros(3, 6)),
            [ -(mu_C*(rmag^2 - 3*x*x))/rmag^5, (3*mu_C*x*y)/rmag^5, (3*mu_C*x*z)/rmag^5, 0, 0, 0, 0, 0, 0, -((lambda_v_mag/p)^(1/(p - 1))*(p*L5^2*L4^2 - L5^2*L5^2 - L6^2*L4^2 - L5^2*L6^2 - L6^2*L5^2 - L6^2*L6^2 - L5^2*L4^2 + p*L5^2*L5^2 + p*L6^2*L4^2 + p*L5^2*L6^2 + p*L6^2*L5^2 + p*L6^2*L6^2 + L4^3*L4 + L4*L5^2*L4 + L4*L6^2*L4))/(1000*mass*(p - 1)*lambda_v_mag^5), -(L4*(1/p)^(1/(p - 1))*lambda_v_mag^2^(1/(2*(p - 1)))*(L5^3 + L4^2*L5 + L5*L6^2 + L5^3 - (p*L5^3) + (L4^2*L5) + (L5*L6^2) - (p*L4^2*L5) - (p*L5*L6^2)))/(1000*mass*(p - 1)*lambda_v_mag^5), -(L4*(1/p)^(1/(p - 1))*lambda_v_mag^2^(1/(2*(p - 1)))*(L6^3 + L4^2*L6 + L5^2*L6 + L6^3 - (p*L6^3) + (L4^2*L6) + (L5^2*L6) - (p*L4^2*L6) - (p*L5^2*L6)))/(1000*mass*(p - 1)*lambda_v_mag^5)]',
            [ (3*mu_C*y*x)/rmag^5, -(mu_C*(rmag^2 - 3*y*y))/rmag^5, (3*mu_C*y*z)/rmag^5, 0, 0, 0, 0, 0, 0, -(L5*(1/p)^(1/(p - 1))*lambda_v_mag^2^(1/(2*(p - 1)))*(L4^3 + L4*L5^2 + L4*L6^2 + L4^3 - (p*L4^3) + (L4*L5^2) + (L4*L6^2) - (p*L4*L5^2) - (p*L4*L6^2)))/(1000*mass*(p - 1)*lambda_v_mag^5), -((lambda_v_mag/p)^(1/(p - 1))*(p*L4^2*L4^2 - L4^2*L5^2 - L4^2*L6^2 - L6^2*L4^2 - L6^2*L5^2 - L6^2*L6^2 - L4^2*L4^2 + p*L4^2*L5^2 + p*L4^2*L6^2 + p*L6^2*L4^2 + p*L6^2*L5^2 + p*L6^2*L6^2 + L5^3*L5 + L4^2*L5*L5 + L5*L6^2*L5))/(1000*mass*(p - 1)*lambda_v_mag^5), -(L5*(1/p)^(1/(p - 1))*lambda_v_mag^2^(1/(2*(p - 1)))*(L6^3 + L4^2*L6 + L5^2*L6 + L6^3 - (p*L6^3) + (L4^2*L6) + (L5^2*L6) - (p*L4^2*L6) - (p*L5^2*L6)))/(1000*mass*(p - 1)*lambda_v_mag^5)]',
            [ (3*mu_C*z*x)/rmag^5, (3*mu_C*z*y)/rmag^5, -(mu_C*(rmag^2 - 3*z*z))/rmag^5, 0, 0, 0, 0, 0, 0, -(L6*(1/p)^(1/(p - 1))*lambda_v_mag^2^(1/(2*(p - 1)))*(L4^3 + L4*L5^2 + L4*L6^2 + L4^3 - (p*L4^3) + (L4*L5^2) + (L4*L6^2) - (p*L4*L5^2) - (p*L4*L6^2)))/(1000*mass*(p - 1)*lambda_v_mag^5), -(L6*(1/p)^(1/(p - 1))*lambda_v_mag^2^(1/(2*(p - 1)))*(L5^3 + L4^2*L5 + L5*L6^2 + L5^3 - (p*L5^3) + (L4^2*L5) + (L5*L6^2) - (p*L4^2*L5) - (p*L5*L6^2)))/(1000*mass*(p - 1)*lambda_v_mag^5), -((lambda_v_mag/p)^(1/(p - 1))*(p*L4^2*L4^2 - L4^2*L5^2 - L5^2*L4^2 - L4^2*L6^2 - L5^2*L5^2 - L5^2*L6^2 - L4^2*L4^2 + p*L4^2*L5^2 + p*L5^2*L4^2 + p*L4^2*L6^2 + p*L5^2*L5^2 + p*L5^2*L6^2 + L6^3*L6 + L4^2*L6*L6 + L5^2*L6*L6))/(1000*mass*(p - 1)*lambda_v_mag^5)]',
            [ (3*mu_C*(- 2*L4*x^3 - 4*L5*x^2*y - 4*L6*x^2*z + 3*L4*x*y^2 + 3*L4*x*z^2 + L5*y^3 + L6*y^2*z + L5*y*z^2 + L6*z^3))/rmag^7, (3*mu_C*(L5*x^3 - 4*L4*x^2*y - 4*L5*x*y^2 - 5*L6*x*y*z + L5*x*z^2 + L4*y^3 + L4*y*z^2))/rmag^7, (3*mu_C*(L6*x^3 - 4*L4*x^2*z + L6*x*y^2 - 5*L5*x*y*z - 4*L6*x*z^2 + L4*y^2*z + L4*z^3))/rmag^7, 0, 0, 0, 0, 0, 0, -(mu_C*(- 2*x^2 + y^2 + z^2))/rmag^5, (3*mu_C*x*y)/rmag^5, (3*mu_C*x*z)/rmag^5]',
            [ (3*mu_C*(L5*x^3 - 4*L4*x^2*y - 4*L5*x*y^2 - 5*L6*x*y*z + L5*x*z^2 + L4*y^3 + L4*y*z^2))/rmag^7, (3*mu_C*(L4*x^3 + 3*L5*x^2*y + L6*x^2*z - 4*L4*x*y^2 + L4*x*z^2 - 2*L5*y^3 - 4*L6*y^2*z + 3*L5*y*z^2 + L6*z^3))/rmag^7, (3*mu_C*(L6*x^2*y + L5*x^2*z - 5*L4*x*y*z + L6*y^3 - 4*L5*y^2*z - 4*L6*y*z^2 + L5*z^3))/rmag^7, 0, 0, 0, 0, 0, 0, (3*mu_C*x*y)/rmag^5, -(mu_C*(x^2 - 2*y^2 + z^2))/rmag^5, (3*mu_C*y*z)/rmag^5]',
            [ (3*mu_C*(L6*x^3 - 4*L4*x^2*z + L6*x*y^2 - 5*L5*x*y*z - 4*L6*x*z^2 + L4*y^2*z + L4*z^3))/rmag^7, (3*mu_C*(L6*x^2*y + L5*x^2*z - 5*L4*x*y*z + L6*y^3 - 4*L5*y^2*z - 4*L6*y*z^2 + L5*z^3))/rmag^7, (3*mu_C*(L4*x^3 + L5*x^2*y + 3*L6*x^2*z + L4*x*y^2 - 4*L4*x*z^2 + L5*y^3 + 3*L6*y^2*z - 4*L5*y*z^2 - 2*L6*z^3))/rmag^7, 0, 0, 0, 0, 0, 0, (3*mu_C*x*z)/rmag^5, (3*mu_C*y*z)/rmag^5, -(mu_C*(x^2 + y^2 - 2*z^2))/rmag^5]',
            hcat(zeros(3, 6), eye(3), zeros(3,3)) )

        if umag > thrustLimit_reduced
            umag = thrustLimit_reduced

            #For the case when thrust is at thrustLimit
            A = vcat(hcat(zeros(3,3), eye(3), zeros(3, 6)),
                [ -(mu_C*(rmag^2 - 3*x*x))/rmag^5, (3*mu_C*x*y)/rmag^5, (3*mu_C*x*z)/rmag^5, 0, 0, 0, 0, 0, 0, -(thrustLimit_reduced*(L5^2 + L6^2))/(1000*mass*lambda_v_mag^3), (L4*L5*thrustLimit_reduced)/(1000*mass*lambda_v_mag^3), (L4*L6*thrustLimit_reduced)/(1000*mass*lambda_v_mag^3)]',
                [ (3*mu_C*y*x)/rmag^5, -(mu_C*(rmag^2 - 3*y*y))/rmag^5, (3*mu_C*y*z)/rmag^5, 0, 0, 0, 0, 0, 0, (L4*L5*thrustLimit_reduced)/(1000*mass*lambda_v_mag^3), -(thrustLimit_reduced*(L4^2 + L6^2))/(1000*mass*lambda_v_mag^3), (L5*L6*thrustLimit_reduced)/(1000*mass*lambda_v_mag^3)]',
                [ (3*mu_C*z*x)/rmag^5, (3*mu_C*z*y)/rmag^5, -(mu_C*(rmag^2 - 3*z*z))/rmag^5, 0, 0, 0, 0, 0, 0, (L4*L6*thrustLimit_reduced)/(1000*mass*lambda_v_mag^3), (L5*L6*thrustLimit_reduced)/(1000*mass*lambda_v_mag^3), -(thrustLimit_reduced*(L4^2 + L5^2))/(1000*mass*lambda_v_mag^3)]',
                [ (3*mu_C*(- 2*L4*x^3 - 4*L5*x^2*y - 4*L6*x^2*z + 3*L4*x*y^2 + 3*L4*x*z^2 + L5*y^3 + L6*y^2*z + L5*y*z^2 + L6*z^3))/rmag^7, (3*mu_C*(L5*x^3 - 4*L4*x^2*y - 4*L5*x*y^2 - 5*L6*x*y*z + L5*x*z^2 + L4*y^3 + L4*y*z^2))/rmag^7, (3*mu_C*(L6*x^3 - 4*L4*x^2*z + L6*x*y^2 - 5*L5*x*y*z - 4*L6*x*z^2 + L4*y^2*z + L4*z^3))/rmag^7, 0, 0, 0, 0, 0, 0, -(mu_C*(- 2*x^2 + y^2 + z^2))/rmag^5, (3*mu_C*x*y)/rmag^5, (3*mu_C*x*z)/rmag^5]',
                [ (3*mu_C*(L5*x^3 - 4*L4*x^2*y - 4*L5*x*y^2 - 5*L6*x*y*z + L5*x*z^2 + L4*y^3 + L4*y*z^2))/rmag^7, (3*mu_C*(L4*x^3 + 3*L5*x^2*y + L6*x^2*z - 4*L4*x*y^2 + L4*x*z^2 - 2*L5*y^3 - 4*L6*y^2*z + 3*L5*y*z^2 + L6*z^3))/rmag^7, (3*mu_C*(L6*x^2*y + L5*x^2*z - 5*L4*x*y*z + L6*y^3 - 4*L5*y^2*z - 4*L6*y*z^2 + L5*z^3))/rmag^7, 0, 0, 0, 0, 0, 0, (3*mu_C*x*y)/rmag^5, -(mu_C*(x^2 - 2*y^2 + z^2))/rmag^5, (3*mu_C*y*z)/rmag^5]',
                [ (3*mu_C*(L6*x^3 - 4*L4*x^2*z + L6*x*y^2 - 5*L5*x*y*z - 4*L6*x*z^2 + L4*y^2*z + L4*z^3))/rmag^7, (3*mu_C*(L6*x^2*y + L5*x^2*z - 5*L4*x*y*z + L6*y^3 - 4*L5*y^2*z - 4*L6*y*z^2 + L5*z^3))/rmag^7, (3*mu_C*(L4*x^3 + L5*x^2*y + 3*L6*x^2*z + L4*x*y^2 - 4*L4*x*z^2 + L5*y^3 + 3*L6*y^2*z - 4*L5*y*z^2 - 2*L6*z^3))/rmag^7, 0, 0, 0, 0, 0, 0, (3*mu_C*x*z)/rmag^5, (3*mu_C*y*z)/rmag^5, -(mu_C*(x^2 + y^2 - 2*z^2))/rmag^5]',
                hcat(zeros(3, 6), eye(3), zeros(3,3)) )
        end
    else
        error("Invalid value of p!")
    end

    #Control vector, scaled properly:
    control_accel = -umag * lambda_v ./ norm(lambda_v) / mass / 1e3 #(kg*m/s^2) -> (km/s^2) acceleration
    if isnan(control_accel[1])
        control_accel = zeros(3)

        #For the case when there is no control (p=1, switching function off)
        A = vcat(hcat(zeros(3,3), eye(3), zeros(3, 6)),
            [-(mu_C*(rmag^2 - 3*x*x))/rmag^5, (3*mu_C*x*y)/rmag^5, (3*mu_C*x*z)/rmag^5, 0, 0, 0, 0, 0, 0, 0, 0, 0]',
            [(3*mu_C*y*x)/rmag^5, -(mu_C*(rmag^2 - 3*y*y))/rmag^5, (3*mu_C*y*z)/rmag^5, 0, 0, 0, 0, 0, 0, 0, 0, 0]',
            [(3*mu_C*z*x)/rmag^5, (3*mu_C*z*y)/rmag^5, -(mu_C*(rmag^2 - 3*z*z))/rmag^5, 0, 0, 0, 0, 0, 0, 0, 0, 0]',
            [(3*mu_C*(- 2*L4*x^3 - 4*L5*x^2*y - 4*L6*x^2*z + 3*L4*x*y^2 + 3*L4*x*z^2 + L5*y^3 + L6*y^2*z + L5*y*z^2 + L6*z^3))/rmag^7, (3*mu_C*(L5*x^3 - 4*L4*x^2*y - 4*L5*x*y^2 - 5*L6*x*y*z + L5*x*z^2 + L4*y^3 + L4*y*z^2))/rmag^7, (3*mu_C*(L6*x^3 - 4*L4*x^2*z + L6*x*y^2 - 5*L5*x*y*z - 4*L6*x*z^2 + L4*y^2*z + L4*z^3))/rmag^7, 0, 0, 0, 0, 0, 0, -(mu_C*(- 2*x^2 + y^2 + z^2))/rmag^5, (3*mu_C*x*y)/rmag^5, (3*mu_C*x*z)/rmag^5]',
            [(3*mu_C*(L5*x^3 - 4*L4*x^2*y - 4*L5*x*y^2 - 5*L6*x*y*z + L5*x*z^2 + L4*y^3 + L4*y*z^2))/rmag^7, (3*mu_C*(L4*x^3 + 3*L5*x^2*y + L6*x^2*z - 4*L4*x*y^2 + L4*x*z^2 - 2*L5*y^3 - 4*L6*y^2*z + 3*L5*y*z^2 + L6*z^3))/rmag^7, (3*mu_C*(L6*x^2*y + L5*x^2*z - 5*L4*x*y*z + L6*y^3 - 4*L5*y^2*z - 4*L6*y*z^2 + L5*z^3))/rmag^7, 0, 0, 0, 0, 0, 0, (3*mu_C*x*y)/rmag^5, -(mu_C*(x^2 - 2*y^2 + z^2))/rmag^5, (3*mu_C*y*z)/rmag^5]',
            [(3*mu_C*(L6*x^3 - 4*L4*x^2*z + L6*x*y^2 - 5*L5*x*y*z - 4*L6*x*z^2 + L4*y^2*z + L4*z^3))/rmag^7, (3*mu_C*(L6*x^2*y + L5*x^2*z - 5*L4*x*y*z + L6*y^3 - 4*L5*y^2*z - 4*L6*y*z^2 + L5*z^3))/rmag^7, (3*mu_C*(L4*x^3 + L5*x^2*y + 3*L6*x^2*z + L4*x*y^2 - 4*L4*x*z^2 + L5*y^3 + 3*L6*y^2*z - 4*L5*y*z^2 - 2*L6*z^3))/rmag^7, 0, 0, 0, 0, 0, 0, (3*mu_C*x*z)/rmag^5, (3*mu_C*y*z)/rmag^5, -(mu_C*(x^2 + y^2 - 2*z^2))/rmag^5]',
            hcat(zeros(3, 6), eye(3), zeros(3,3)) )

    end


    #Phi-dot:
    ϕ̇ = A * ϕ

    #Output in place: (scaled)
    dstate_costate[1:3] = state_costate_unscaled[4:6] ./ scale_vec[1:3]
    dstate_costate[4:6] = (-mu_C/norm(state_costate_unscaled[1:3])^3 * state_costate_unscaled[1:3] + control_accel) ./ scale_vec[4:6]

    dstate_costate[7] = ((mu_C*(2*L4*x.^2 + 3*L5*x*y + 3*L6*x*z - L4*y.^2 - L4*z.^2))./(rmag)^5) / scale_vec[7]
    dstate_costate[8] = ((mu_C*(- L5*x.^2 + 3*L4*x*y + 2*L5*y.^2 + 3*L6*y*z - L5*z.^2))./(rmag)^5) / scale_vec[8]
    dstate_costate[9] = ((mu_C*(- L6*x.^2 + 3*L4*x*z - L6*y.^2 + 3*L5*y*z + 2*L6*z.^2))./(rmag)^5) / scale_vec[9]
    dstate_costate[10:12] = state_costate_unscaled[7:9] ./ scale_vec[10:12]

    dstate_costate[13:end] = ϕ̇[:] ./ scale_mat_vec #re-scale
end
