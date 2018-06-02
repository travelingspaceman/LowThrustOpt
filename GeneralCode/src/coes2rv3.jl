#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#COES to r,v
#all angles in degrees
#
#In:
#   a = semimajor axis
#   e = eccentricity
#   inc = inclination
#   RAAN = right ascension of ascending node
#   arg_per = argument of perigee
#   theta = true anomaly
#   mu = gravitational parameter
#Out:
#   r = position vector in km
#   v = velocity vector in km/s

function coes2rv3(sma, ecc, inc, RAAN, AOP, TA, mu_C)

  h = sqrt(sma*mu_C*(1-ecc^2));

  r_x = h^2/mu_C/(1+ecc*cosd(TA))*[cosd(TA); sind(TA); 0]; #perifocal frame
  v_x = mu_C/h*[-sind(TA); ecc+cosd(TA); 0]; #perifocal frame

  Q_xX = geo_peri(RAAN, inc, AOP)

  R = Q_xX * r_x
  V = Q_xX * v_x

  #Outputs:
  (R, V)
end



function geo_peri(RAAN,inc,arg_per)
  #Coordinate transformation between geocentric equatorial and perifocal
  #frames
  # X = perifocal
  # x = geocentric equatorial

  #row 1
  Q_Xx_1 = [-sind(RAAN)*cosd(inc)*sind(arg_per)+cosd(RAAN)*cosd(arg_per);
      cosd(RAAN)*cosd(inc)*sind(arg_per)+sind(RAAN)*cosd(arg_per);  sind(inc)*sind(arg_per)]

  #row 2
  Q_Xx_2 = [ -sind(RAAN)*cosd(inc)*cosd(arg_per)-cosd(RAAN)*sind(arg_per);
      cosd(RAAN)*cosd(inc)*cosd(arg_per)-sind(RAAN)*sind(arg_per);  sind(inc)*cosd(arg_per)]

  #row 3
  Q_Xx_3 = [sind(RAAN)*sind(inc);  -cosd(RAAN)*sind(inc);  cosd(inc)]

  #output:
  Q_Xx = hcat(Q_Xx_1, Q_Xx_2, Q_Xx_3)
end
