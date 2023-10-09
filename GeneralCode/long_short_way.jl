#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Determines long way or short way for Lambert's problem solutions. Chooses
#whichever solution gives a prograde orbit
#
#Inputs:
#   r_pd = R-vector of departing planet (heliocentric coord.)
#   r_pa = R-vector of arriving planet (heliocentric coord.)
#Outputs:
#   lw = 0 for short way
#      = 1 for long way

function long_short_way(r_pd,r_pa)

  angle1 = atan2(r_pd[2],r_pd[1]);
  angle2 = atan2(r_pa[2],r_pa[1]);

  if angle1<0 #put angle1 between 0->2*pi
      angle1 = angle1+2*pi;
  end

  if angle2<0 #put angle2 between 0->2*pi
      angle2 = angle2+2*pi;
  end
  temp = angle2-angle1;

  #initialize:
  lw = false

  if temp > 0
      if temp < pi
          lw=false; #short way
      elseif temp > pi
          lw=true; #long way
      end
  elseif temp < 0
      if temp < -pi
          lw=false; #short way
      elseif temp > -pi
          lw=true; #long way
      end
  end

  #Output:
  lw
end
