#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Linear interpolation, 1-dimensional
#
#Make sure that xi is within the range of x -- no checks are included here

function LinInterp(x,y,xi)
  #Initialize values that we find in the for loop:
  #if xi is at the beginning of x:
  y0 = y[(1)];
  y1 = y[(2)];
  x0 = x[(1)];
  x1 = x[(2)];

  for i = 1:length(x)
    if x[i] > xi
      x0 = x[i-1]
      x1 = x[i]

      y0 = y[i-1]
      y1 = y[i]

      break
    end

    #if xi is at the end of x:
    x0 = x[(end-1)];
    x1 = x[(end)];
    y0 = y[(end-1)];
    y1 = y[(end)];
  end

  yi = y0 + (xi-x0)*(y1-y0)/(x1-x0)
end
