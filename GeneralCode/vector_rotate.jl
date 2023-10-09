#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Rotate any vector about any other vector
#rotation math from:
# http://steve.hollasch.net/cgindex/math/rotvec.html
#
#Inputs:
#   V_minus  = 3-element vector before rotation
#   rot_axis = 3-element rotation axis
#   alpha    = angle to rotate vector by CCW (radians)
#Outputs:
#   V_plus   = vector after rotation [1x3]

function vector_rotate(V_minus, rot_axis, alpha)

  L = [0 rot_axis[3] -rot_axis[2];
      -rot_axis[3] 0 rot_axis[1];
      rot_axis[2] -rot_axis[1] 0]; #some cool matrix

  d = norm(rot_axis); #magnitude of vector to rotate about

  V_plus = (eye(3) + sin(alpha)/d*L + ((1-cos(alpha))/(d^2)*(L*L)))' * V_minus

end
