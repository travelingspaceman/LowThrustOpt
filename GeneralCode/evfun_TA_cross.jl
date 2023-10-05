#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Events functions -- stop at true anomaly

function evfun_TA_cross(t, state, varargin...)
  #detect true anomaly crossing

  mu_C = varargin[3]
  time_direction = varargin[4]
  TA_stop = varargin[5]

  #First, compute the orbital elements
  (sma,ecc,inc,RAAN,AOP,TA,h,P,E,e_vec,h_vec) = rv2coes(state[1:3], state[4:6], 1, mu_C) #deg


  #set direction to forward
  dirn = 1 #from (-) to (+)

  dirn *= time_direction #modify event direction for backwards propagation

  #events triggered when evf == 0
  evf = TA_stop - TA

  #outputs:
  evf, dirn
end
