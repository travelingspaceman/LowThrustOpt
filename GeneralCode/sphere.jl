#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory

function sphere(N)

  #create data points for a sphere surface:
  u = LinRange(0, 2π, N)
  v = LinRange(0, π, N)
  x = cos.(u) * sin.(v)'
  y = sin.(u) * sin.(v)'
  z = repeat(cos.(v)',outer=[N, 1])

  #outputs:
  (x,y,z)

end
