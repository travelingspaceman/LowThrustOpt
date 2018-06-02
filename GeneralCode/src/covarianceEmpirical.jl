#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory
#
#Compute the empirical covariance matrix
#
#Inputs:
#   x       = (m x N) array of samples. N samples, where each sample is m elements
#Outputs:
#   P       = (m x N) empirical covariance matrix
#   xbar    = (m) average vector

function covarianceEmpirical(x)
    N = size(x, 2)
    xbar = mean(x,2)
    P = 1/N * (x - repmat(xbar, 1, N)) * (x - repmat(xbar, 1, N))'

    (P, xbar)
end
