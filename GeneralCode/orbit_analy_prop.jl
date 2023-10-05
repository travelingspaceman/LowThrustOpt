#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Two-body orbit analytical propagation (elliptical orbit only).
#
#Inputs:
#   R0      = Initial position vector [3x1] (km)
#   V0      = Initial velocity vector [3x1] (km/s)
#   dt      = Amount of time to propagate [N] (sec)
#   mu_C    = gravitational parameter of central body (km3/s2)
#Outpus:
#   Rf      = Final position vector [3xN] (km)
#   Vf      = Final velocity vector [3xN] (km)

function orbit_analy_prop(R0, V0, dt, mu_C)
    ## initialize outputs
    Rf = NaN*ones(3,length(dt));
    Vf = NaN*ones(3,length(dt));
    TA_deg = NaN*ones(1,length(dt));

    ## First, find orbital elements of initial state:
    (sma,ecc,inc,RAAN,arg_per,TA0,h,Period,E,e_vec,h_vec) = rv2coes(R0, V0, 2, mu_C)
    if isnan(RAAN)
        RAAN = 0;
    end
    if isnan(arg_per)
        arg_per = 0;
    end
    if ecc >= 1
        return (NaN, NaN, NaN) #will return NaN's
    end
    n = 2*pi / Period;

    #compute Eccentric anomaly:
    E0 = acos.( (ecc + cos.(TA0))/(1 + ecc*cos.(TA0)) ); #rad
    if TA0 > pi
        E0 = 2*pi - E0;
    end

    #mean anomaly:
    M0 = E0 - ecc*sin.(E0); #rad

    #time past periapse:
    tp0 = M0 / n;

    ## find orbital elements  of final state:
    Mf = n* (tp0 + dt);
    Mf = mod.(Mf, 2*pi);

    # initial guess for eccentric anomaly:
    if length(dt) == 1
        if Mf < pi
            Ef = Mf + ecc/2
        else
            Ef = Mf - ecc/2
        end
    else
        Ef = zeros(size(Mf))
        for ind = 1:length(dt)
            if Mf[ind] < pi
                Ef[ind] = Mf[ind] + ecc/2
            else
                Ef[ind] = Mf[ind] - ecc/2
            end
        end
    end

    # iterate using Newton's Method:
    tol = 1e-12;
    err = 1;
    num_iter = 0;
    while err > tol
        f = Ef - ecc*sin.(Ef) - Mf;
        f1 = 1-ecc*cos.(Ef); #derivative of f (wrt E)

        Ef_new = Ef - f ./ f1;
        err = norm(Ef_new - Ef, Inf);

        Ef = Ef_new;

        num_iter = num_iter+1;
        if num_iter > 20
            return #will return NaN's
        end
    end


    #use 'atan' instead of 'atan2' to be friendly for Dual Numbers:
    # TAf = 2*atan2.(tan.(Ef/2), sqrt((1-ecc)/(1+ecc)))
    TAf = atan.( tan.(Ef/2) / sqrt((1-ecc)/(1+ecc)) )
    #quadrant check:
    if sqrt((1-ecc)/(1+ecc)) < 0
        TAf += pi
    end
    TAf *= 2

    ## find final state [Rf, Vf]
    #rotation matrices:
    Q_Xx = vcat([-sin.(RAAN)*cos.(inc)*sin.(arg_per)+cos.(RAAN)*cos.(arg_per),
        cos.(RAAN)*cos.(inc)*sin.(arg_per)+sin.(RAAN)*cos.(arg_per), sin.(inc)*sin.(arg_per)]',
        [-sin.(RAAN)*cos.(inc)*cos.(arg_per)-cos.(RAAN)*sin.(arg_per),
        cos.(RAAN)*cos.(inc)*cos.(arg_per)-sin.(RAAN)*sin.(arg_per), sin.(inc)*cos.(arg_per)]',
        [sin.(RAAN)*sin.(inc), -cos.(RAAN)*sin.(inc), cos.(inc)]')
    Q_xX = Q_Xx';

    #dual numbers friendly:
    Rf = eltype(R0).(ones(3,length(dt)))
    Vf = eltype(R0).(ones(3,length(dt)))

    for ind = 1:length(dt)
        #define position / velocity vectors in perifocal frame
        r_x = h^2 ./ mu_C ./ (1+ecc*cos.(TAf[ind])) .* [cos.(TAf[ind]), sin.(TAf[ind]), 0.]; #perifocal frame
        v_x = mu_C/h * [-sin.(TAf[ind]), ecc+cos.(TAf[ind]), 0.]; #perifocal frame

        #rotate into inertial frame:
        Rf[:,ind] = Q_xX * r_x;
        Vf[:,ind] = Q_xX * v_x;
    end

    TA_deg = TAf * 180/pi;

    #Outputs:
    (Rf, Vf, TA_deg)
end
