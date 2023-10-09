#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Converts from Cartesian coordinates (position, velocity) to Modified
#Equinoctial Elements

function rv2mee(R, V, mu_C)

  # convert eci state vector to modified equinoctial elements

  # input

  #  mu_C   = gravitational constant (km**3/sec**2)
  #  reci = eci position vector (kilometers)
  #  veci = eci velocity vector (kilometers/second)

  # output

  #  mee(1) = semiparameter (kilometers)
  #  mee(2) = f equinoctial element
  #  mee(3) = g equinoctial element
  #  mee(4) = h equinoctial element
  #  mee(5) = k equinoctial element
  #  mee(6) = true longitude (radians)

  # Orbital Mechanics with Matlab

  ###############################

  radius = norm(R);

  hv = cross(R, V);

  hmag = norm(hv);

  pmee = hmag^2 / mu_C;

  rdotv = dot(R, V);

  rzerod = rdotv / radius;

  eccen = cross(V, hv);

  uhat = R / radius;
  vhat = (radius * V - rzerod * R) / hmag;

  eccen = eccen / mu_C - uhat;

  # unit angular momentum vector
  hhat = hv / norm(hv);

  # compute kmee and hmee
  denom = 1.0 + hhat[3];
  kmee = hhat[1] / denom;
  hmee = -hhat[2] / denom;

  # construct unit vectors in the equinoctial frame
  fhat = eltype(R).(zeros(3))
  fhat[1] = 1.0 - kmee^2 + hmee^2;
  fhat[2] = 2.0 * kmee * hmee;
  fhat[3] = -2.0 * kmee;

  ghat = eltype(R).(zeros(3))
  ghat[1] = fhat[2];
  ghat[2] = 1.0 + kmee^2 - hmee^2;
  ghat[3] = 2.0 * hmee;

  ssqrd = 1.0 + kmee^2 + hmee^2;

  # normalize
  fhat = fhat / ssqrd;
  ghat = ghat / ssqrd;

  # compute fmee and gmee
  fmee = dot(eccen, fhat);
  gmee = dot(eccen, ghat);

  # compute true longitude
  cosl = uhat[1] + vhat[2];
  sinl = uhat[2] - vhat[1];
  lmee = atan2(sinl, cosl);

  p = pmee;
  f = fmee;
  g = gmee;
  h = hmee;
  k = kmee;
  L = lmee;

  MEE = [p; f; g; h; k; L] #all in one vector for convenience

  #outputs:
  (MEE, p, f, g, h, k, L)
end



function rv2mee_r(R, V, mu_C, direction)
  #Converts from Cartesian coordinates (position, velocity) to Modified
  #Equinoctial Elements, for retrograde element set. This moves the
  #singularity from inc=180deg to inc=0deg.

  #convert to COE's:
  (sma,ecc,inc,RAAN,AOP,TA,h,P,E,e_vec,h_vec) = rv2coes(R,V,2,mu_C) #km, rad

  if isnan(RAAN)
    RAAN = 0.0
  end
  if isnan(AOP)
    AOP = 0.0
  end
  if isnan(TA)
    TA = 0.0
  end

  p = h^2 / mu_C; #semilatus rectum

  f = ecc * cos(AOP + direction * RAAN);

  g = ecc * sin(AOP + direction * RAAN);

  if direction == 1
      h = tan(inc/2) * cos(RAAN);
  elseif direction == -1
      h = atan(inc/2) * cos(RAAN);
  end

  if direction == 1
      k = tan(inc/2) * sin(RAAN);
  elseif direction == -1
      k = atan(inc/2) * sin(RAAN);
  end

  L = AOP + direction*RAAN + TA;

  MEE = [p; f; g; h; k; L] #all in one vector for convenience

  #invalidate results if really close to singularities:
  if (abs(inc - pi) < 1e-3) & (direction == 1)
    MEE = NaN*MEE
  elseif (abs(inc) < 1e-3) & (direction == -1)
    MEE = NaN*MEE
  end

  #Outputs:
  (MEE, p, f, g, h, k, L)

end
