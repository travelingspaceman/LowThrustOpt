#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Contains funcsions 'synodic2inertial' and 'inertial2synodic' to transfer
#between synodic (rotating) and inertial frames in the CRTBP.

function synodic2inertial(R_synodic, V_synodic, t_sec, DU, TU, R_central)
  #Inputs:
  # R_synodic -- in the synodic reference frame, dimensionless units
  # V_synodic
  # t_sec -- seconds since t0
  #
  #Outputs:
  # R_sc_inr_central_km -- in inertial frame, relative to "central" body, km
  # V_sc_inr_central_km

  #rotation of synodic frame compared to inertial frame:
  theta0 = 0.
  thetadot = 1 / TU #orbit rotation rate of Earth-Moon system (rad/sec)

  omega = [0; 0; thetadot]

  #Convert to dimensional units:
  #position:
  R_sc_rot_bary_km = R_synodic * DU; # (DU) -> (km)
  #velocity:
  V_sc_rot_bary_km = V_synodic * DU/TU; # (DU/TU) -> (km/s)

  #angle between frames:
  theta = theta0 + thetadot * t_sec; #rad

  #unit vectors defining rotating coordinate frame
  r = [cos(theta); sin(theta); 0] #direction, central -> 3rd body
  x = r ./ norm(r)
  z = omega ./ norm(omega)
  y = cross(z, x)

  #direction cosine matrix:
  C = hcat(x, y, z) #rotating frame --> inertial frame
  Cdot = skewSymmetric(omega) * C
  #Helpful definition:
  #   Cdot*R === cross(omega, C*R)

  #position of spacecraft in rotating frame, relative to central body
  R_sc_rot_central_km = R_sc_rot_bary_km - R_central;

  #Rotate r,v from rotating frame to inertial frame:
  R_sc_inr_central_km = C * R_sc_rot_central_km;
  V_sc_inr_central_km = C * V_sc_rot_bary_km + Cdot * R_sc_rot_central_km;

  #Outputs:
  (R_sc_inr_central_km, V_sc_inr_central_km, C)
end


function inertial2synodic(R_sc_inr_central_km, V_sc_inr_central_km, t_sec, DU, TU, R_central)
  #Inputs:
  # R_synodic -- in the synodic reference frame, dimensionless units
  # V_synodic
  # t_sec -- seconds since t0
  #
  #Outputs:
  # R_sc_inr_central_km -- in inertial frame, relative to "central" body, km
  # V_sc_inr_central_km

  #rotation of synodic frame compared to inertial frame:
  theta0 = 0.
  thetadot = 1 / TU #orbit rotation rate of Earth-Moon system (rad/sec)

  omega = [0; 0; thetadot]

  #angle between frames:
  theta = theta0 + thetadot * t_sec; #rad

  #unit vectors defining rotating coordinate frame
  r = [cos(theta); sin(theta); 0] #direction, central -> 3rd body
  x = r ./ norm(r)
  z = omega ./ norm(omega)
  y = cross(z, x)

  #direction cosine matrix:
  C = hcat(x, y, z) #rotating frame --> inertial frame
  Cdot = skewSymmetric(omega) * C
  #Helpful definition:
  #   Cdot*R === cross(omega, C*R)


  #Rotate from inertial to rotating frame:
  R_sc_rot_central_km = C' * R_sc_inr_central_km
  V_sc_rot_bary_km = C' * (V_sc_inr_central_km - Cdot * R_sc_rot_central_km)

  #position of spacecraft in rotating frame, relative to barycenter
  R_sc_rot_bary_km = R_sc_rot_central_km + R_central

  #Convert to dimensionless units:
  #position:
  R_synodic = R_sc_rot_bary_km / DU # (km) -> (DU)
  #velocity:
  V_synodic = V_sc_rot_bary_km / (DU/TU) # (km/s) -> (DU/TU)


  #Outputs:
  (R_synodic, V_synodic)
end
