#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Derivatives function for 2-body orbit propagation

function orbit_num_prop_deriv(t, state, varargin...)
  #for use with "ode78()"

  if length(varargin) < 1
    #default to Earth orbit
    mu_C = 398600.4418;
  else
    mu_C = varargin[1]
  end

  mag_R = norm(state[1:3]);

  dstate = zeros(6)
  dstate[1:3] = state[4:6]
  dstate[4] = -mu_C*state[1]./mag_R.^3
  dstate[5] = -mu_C*state[2]./mag_R.^3
  dstate[6] = -mu_C*state[3]./mag_R.^3
end


#For new version of "DifferentialEquations.jl"
function orbit_num_prop_deriv!(dstate, state, params, t)
  #pull out parameters:
  mu_C = params

  mag_R = norm(state[1:3]);

  #update dstate in place:
  dstate[1:3] = state[4:6]
  dstate[4] = -mu_C*state[1]./mag_R.^3
  dstate[5] = -mu_C*state[2]./mag_R.^3
  dstate[6] = -mu_C*state[3]./mag_R.^3
end
