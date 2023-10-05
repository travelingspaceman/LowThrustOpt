#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Two-body propagator with thrust as given


#for use with new version of DifferentialEquations.jl
function TwoBody_prop_EP_deriv!(dstate, state, params, t::Float64)

  #pull out parameters:
  (Isp, control, mu_C, mass0, time_direction, J2_on) = params

  #pull out from state vector:
  x = state[1];
  y = state[2];
  z = state[3];
  if length(state) == 7
      m = state[7];
  else
      m = mass0
  end

  #thrust acceleration magnitude (km/s^2):
  T_mag = norm(control)/m/1e3; # (kg*m/s^2) -> (km/s^2) acceleration

  #thrust acceleration vector:
  if length(control) == 1
    #vector not specified. Thrust in velocity direction
    control = state[4:6] / norm(state[4:6]) * control
  end

  if norm(control) == 0
    T = control
  else
    T = control ./ norm(control) * T_mag
  end


  #calculate accelerations:
  r = norm(state[1:3]);

  if J2_on == true
    #Add J2 perturbation:
    J2 = 1.082626925638815e-3;
    R_planet = 6378.
    xddot = -mu_C*x/r^3*(1-3*J2*R_planet^2/(2*r^2)*(5*z^2/r^2 - 1)) + T[1]; #km/s2
    yddot = -mu_C*y/r^3*(1-3*J2*R_planet^2/(2*r^2)*(5*z^2/r^2 - 1)) + T[2];
    zddot = -mu_C*z/r^3*(1-3*J2*R_planet^2/(2*r^2)*(5*z^2/r^2 - 3)) + T[3];
  else
    xddot = -mu_C*x./r.^3 + T[1];
    yddot = -mu_C*y./r.^3 + T[2];
    zddot = -mu_C*z./r.^3 + T[3];
  end

  #Modify dstate in place:
  dstate[1:3] = state[4:6]
  dstate[4] = xddot
  dstate[5] = yddot
  dstate[6] = zddot
  if length(state) == 7 #mass flow rate, if mass is used
    g0 = 9.80665; #m/s2
    dstate[7] = - time_direction * norm(control)/(Isp*g0); #mass flow rate (kg/s)
  end

end





#for new version of DifferentialEquations
function TwoBody_prop_EP_precomputedControl_deriv!(dstate, state, params, t::Float64)
  #for use with DifferentialEquations.jl
  #
  #Function to propagate an orbit with a pre-computed thrust history.

  #pull out parameters:
  (Isp, mu_C, mass0, time_direction, J2_on, u_itp) = params

  #pull out from state vector:
  x = state[1];
  y = state[2];
  z = state[3];
  if length(state) == 7
      m = state[7];
  else
      m = mass0
  end

  #u_itp is an interpolation object for the control. It is prepared like this:
  # using Interpolations
  # u_itp = scale(interpolate(u_all', (BSpline(Cubic(Natural())), NoInterp()), OnGrid()), t_sec, 1:3)
  # (from https://github.com/JuliaMath/Interpolations.jl)
  #where t_sec has to be an evenly-spaced time grid.

  #interpolate control from given history:
  control = [u_itp[t, 1], u_itp[t, 2], u_itp[t, 3]]

  #u_hist is given in Newtons. Thrust acceleration magnitude (km/s^2):
  T_mag = norm(control)/m/1e3; # (kg*m/s^2) -> (km/s^2) acceleration

  #thrust acceleration vector:
  if norm(control) == 0
    T = control
  else
    T = control ./ norm(control) * T_mag
  end

  #calculate accelerations:
  r = norm(state[1:3]);

  if J2_on == true
    #Add J2 perturbation:
    J2 = 1.082626925638815e-3;
    R_planet = 6378.
    xddot = -mu_C*x/r^3*(1-3*J2*R_planet^2/(2*r^2)*(5*z^2/r^2 - 1)) + T[1]; #km/s2
    yddot = -mu_C*y/r^3*(1-3*J2*R_planet^2/(2*r^2)*(5*z^2/r^2 - 1)) + T[2];
    zddot = -mu_C*z/r^3*(1-3*J2*R_planet^2/(2*r^2)*(5*z^2/r^2 - 3)) + T[3];
  else
    xddot = -mu_C*x./r.^3 + T[1];
    yddot = -mu_C*y./r.^3 + T[2];
    zddot = -mu_C*z./r.^3 + T[3];
  end

  #Modify dstate in place:
  dstate[1:3] = state[4:6]
  dstate[4] = xddot
  dstate[5] = yddot
  dstate[6] = zddot
  if length(state) == 7 #mass flow rate, if mass is used
    g0 = 9.80665; #m/s2
    dstate[7] = - time_direction * norm(control)/(Isp*g0); #mass flow rate (kg/s)
  end

end


#for new version of DifferentialEquations
function TwoBody_prop_EP_NNControl_deriv!(dstate, state, params, t::Float64)
  #for use with DifferentialEquations.jl
  #
  #Function to propagate an orbit with a neural-network updated control

  #pull out parameters:
  (Isp, mu_C, AU, mass0, time_direction, J2_on, λv_itp, PowerSystem, thrustFactor, p, ρ) = params

  #pull out from state vector:
  x = state[1];
  y = state[2];
  z = state[3];
  if length(state) == 7
      m = state[7];
  else
      m = mass0
  end

  #λv_itp is an interpolation object for the control. It is prepared like this:
  # using Interpolations
  # λv_itp = scale(interpolate(λv_all', (BSpline(Cubic(Natural())), NoInterp()), OnGrid()), t_sec, 1:3)
  # (from https://github.com/JuliaMath/Interpolations.jl)
  #where t_sec has to be an evenly-spaced time grid.

  #interpolate control from given history:
  λv = [λv_itp[t, 1], λv_itp[t, 2], λv_itp[t, 3]]

  #power available at time t for thrusting:
  thrustLimit_reduced = PowerSystem( t / 86400 / 365.25, norm(state[1:3] / AU), thrustFactor )


  #Thrust magnitude from control law:
  if p == 1 #minimum mass problem
      #Junkins' approach: smooth the discontinuity
      g = norm(λv) - 1
      umag = 1/2 * (thrustLimit_reduced + thrustLimit_reduced*tanh(g/ρ) )

  elseif p > 1 #p is in the range (1, 2]
      umag = (1/p * norm(λv))^(1 / (p-1));

      if umag > thrustLimit_reduced
          umag = thrustLimit_reduced
      end
  else
      error("Invalid value of p!")
  end

  control = -λv ./ norm(λv) .* umag #3-vector



  #u_hist is given in Newtons. Thrust acceleration magnitude (km/s^2):
  T_mag = norm(control)/m/1e3; # (kg*m/s^2) -> (km/s^2) acceleration

  #thrust acceleration vector:
  if norm(control) == 0
    T = control
  else
    T = control ./ norm(control) * T_mag
  end

  #calculate accelerations:
  r = norm(state[1:3]);

  if J2_on == true
    #Add J2 perturbation:
    J2 = 1.082626925638815e-3;
    R_planet = 6378.
    xddot = -mu_C*x/r^3*(1-3*J2*R_planet^2/(2*r^2)*(5*z^2/r^2 - 1)) + T[1]; #km/s2
    yddot = -mu_C*y/r^3*(1-3*J2*R_planet^2/(2*r^2)*(5*z^2/r^2 - 1)) + T[2];
    zddot = -mu_C*z/r^3*(1-3*J2*R_planet^2/(2*r^2)*(5*z^2/r^2 - 3)) + T[3];
  else
    xddot = -mu_C*x./r.^3 + T[1];
    yddot = -mu_C*y./r.^3 + T[2];
    zddot = -mu_C*z./r.^3 + T[3];
  end

  #Modify dstate in place:
  dstate[1:3] = state[4:6]
  dstate[4] = xddot
  dstate[5] = yddot
  dstate[6] = zddot
  if length(state) == 7 #mass flow rate, if mass is used
    g0 = 9.80665; #m/s2
    dstate[7] = - time_direction * norm(control)/(Isp*g0); #mass flow rate (kg/s)
  end

end
