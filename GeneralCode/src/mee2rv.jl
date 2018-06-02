#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Converts from Modified Equinoctial Elements to Cartesian coordinates
#
#Can do batch conversion by inputting MEE's as row vectors of length N
#Outputs will then be 3-by-N vectors R and V

#ONLY WORKS FOR PROGRADE ELEMENTS!!!!

function mee2rv(p, f, g, h, k, L, mu_C)

  alpha2 = h.^2 - k.^2; #alpha^2
  s2 = 1 + h.^2 + k.^2;
  w = 1 + f.*cos(L) + g.*sin(L);
  r = p./w;

  R = [ r./s2 .* (cos(L) + alpha2.*cos(L) + 2.*h.*k.*sin(L));
      r./s2 .* (sin(L) - alpha2.*sin(L) + 2.*h.*k.*cos(L));
      2.*r./s2 .* (h.*sin(L) - k.*cos(L))];

  V = [ -1./s2 .* sqrt(mu_C./p) .* (sin(L) + alpha2.*sin(L) - 2.*h.*k.*cos(L) + g - 2.*f.*h.*k + alpha2.*g);
      -1./s2 .* sqrt(mu_C./p) .* (-cos(L) + alpha2.*cos(L) + 2.*h.*k.*sin(L) - f + 2.*g.*h.*k + alpha2.*f);
      2./s2 .* sqrt(mu_C./p) .* (h.*cos(L) + k.*sin(L) + f.*h + g.*k)];

  #outputs
  (R,V)
end
