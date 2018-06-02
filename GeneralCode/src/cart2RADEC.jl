#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Functions to convert between RA/DEC/r and x/y/z

#Convert a Cartesian vector R to Right Ascension and Declination (radians)
function cart2RaDec(R)
    RA = atan2(R[2], R[1]) #rad
    r = norm(R)
    DEC = asin(R[3] / r) #rad

    #outputs:
    [RA, DEC, r]
end

#Convert Right Ascension and Declination to Cartesian vector
function RaDec2cart(RA, DEC, r)
    rxy = r * cos(DEC)
    z = r * sin(DEC)

    x = rxy * cos(RA)
    y = rxy * sin(RA)

    #Outputs:
    [x, y, z]
end
