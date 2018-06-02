#Copyright Nathan Parrish 2018
#University of Colorado, Boulder
#Colorado Center for Astrodynamics Research
#Celestial and Spaceflight Mechanics Laboratory

function norm_many(V)
    #Input: a [m x N] matrix
    #Return: the 2-norm of each column

    vv2 = V[1,:].^2
    for ind = 2:size(V, 1)
        vv2 += V[ind,:].^2
    end

    vv = sqrt.(vv2)
end
