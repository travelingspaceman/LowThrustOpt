#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#This function has inputs of vector R [km] and vector V [km/s] of an
#orbiting body. It outputs the Classical Orbital Elements: semi-major
#axis a, eccentricity e, inclination i, longitude of ascending node OM,
#argument of perigee arg_per, true anomaly nu
#
#only works for ELLIPTICAL ORBITS

function rv2coes(R,V,deg_rad,mu_C)
  # initialize some stuff
  RAAN = NaN;
  AOP = NaN;

  mag_R = norm(R);
  mag_V = norm(V);
  h_vec = cross(R,V); #specific angular momentum (perpendicular to orbit plane
  h = norm(h_vec);

  ## semi-major axis
  E = (mag_V^2)/2 - mu_C/mag_R; #specific mechanical energy
  sma = -mu_C/(2*E);

  ## eccentricity
  e_vec = 1/mu_C*((mag_V^2 - mu_C/mag_R)*R - dot(R,V)*V);
  ecc = norm(e_vec);

  ## inclination
  inc = acos(dot([0;0;1],h_vec)/norm(h_vec));
  if inc >pi
      inc = inc-pi;
  end

  ## longitude of ascending node
  n = cross([0;0;1],h_vec);   #ascending node vector
  RAAN = acos(dot([1;0;0],n)/norm(n));
  if RAAN < pi
      if n[2] < 0
          RAAN = 2*pi - RAAN;
      end
  end
  if RAAN > pi
      if n[2] > 0
          RAAN = 2*pi - RAAN;
      end
  end

  ## argument of perigee
  AOP = acos(dot(n,e_vec)/(norm(n)*ecc));
  if (dot(n,e_vec)/(norm(n)*ecc)) > 1
      AOP = 0;
  end
  if AOP < pi
      if e_vec[3] < 0
          AOP = 2*pi - AOP;
      end
  end
  if AOP > pi
      if e_vec[3] > 0
          AOP = 2*pi - AOP;
      end
  end

  ## true anomaly
  temp = dot(e_vec,R)/(ecc*mag_R);
  if temp > 1  #sometimes have rounding errors where temp = 1 + 1e-16
      temp = 1;
  end
  if temp < -1.
    temp = -1.;
  end
  # println(temp)
  TA = acos(temp);
  if TA < pi
      if dot(R,V) < 0
          TA = 2*pi - TA;
      end
  end
  if TA > pi
      if dot(R,V) > 0
          TA = 2*pi - TA;
      end
  end
  TA = real(TA); #sometimes get imaginary TA because of numerical issues

  ## Special cases
  if (inc == 0) && (ecc != 0) #equatorial, non-circular
      RAAN = 0;

      w_true = atan2(e_vec[(2)], e_vec[(1)]);
      if h_vec[(3)] < 0
          w_true = 2*pi - w_true;
      end

      AOP = w_true;
  elseif (inc == pi) && (ecc != 0)
      RAAN = 0;

      w_true = atan2(e_vec[(2)], e_vec[(1)]);
      if h_vec[(3)] < 0
          w_true = 2*pi - w_true;
      end

      AOP = w_true;

  elseif (inc != 0) && (ecc == 0) #inclined, circular
      AOP = 0;

      u = acos(dot(n,R)/(norm(R)*norm(n))); #argument of latitude (arg_per + TA)
      if R[(3)] < 0
          u = 2*pi - u;
      end

      TA = u;
  elseif (inc == 0) && (ecc == 0) #equatorial, circular
      RAAN = 0;
      AOP = 0;

      lambda = acos(R[(1)]/norm(R)); #true longitude
      if R[(2)] < 0
          lambda = 2*pi - lambda;
      end

      TA = lambda;
  end

  ## orbit period
  P = 2*pi*sqrt((sma^3)/mu_C);

  ## convert radians to degrees, if inputs specify
  # if deg_rad == 'd'
  if deg_rad == 1 #degrees
      inc = inc*180/pi;
      RAAN = RAAN*180/pi;
      AOP = AOP*180/pi;
      TA = TA*180/pi;
  # elseif deg_rad == 'r'
  elseif deg_rad == 2 #radians
  #do nothing to stay in radians
  else
      inc = inc*180/pi;
      RAAN = RAAN*180/pi;
      AOP = AOP*180/pi;
      TA = TA*180/pi;
  end

  (sma,ecc,inc,RAAN,AOP,TA,h,P,E,e_vec,h_vec)

end
