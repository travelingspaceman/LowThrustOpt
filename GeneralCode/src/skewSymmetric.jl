#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Converts a [3x1] vector to its [3x3] equivalent skew-symmetric matrix.
# Q_tilde = cross(q,_)
#This is the matrix math equivalent of a cross product

function skewSymmetric(q)

  Q_tilde =
      [  0.  -q[3]    q[2];
       q[3]     0.   -q[1];
      -q[2]   q[1]      0.]

end
