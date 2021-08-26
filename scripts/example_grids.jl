#=
Created on Thursday 30 April 2020
Last update: Friday 1 May 2020

@author: Michiel Stock
michielfmstock@gmail.com

Script for bachelor students, explaining the grid init and contour plot.
=#

# initializing a bacteriagrid
# ---------------------------

using PhageSim

# generate grid of size 100 by 100, containing 200 bacteria over three species
mybactgrid = initbactgrid(100, 100, nbacteria=300, nspecies=3)

#check, density of bacteria, i.e. numb of bact / number of spaces
density(mybactgrid)  # 0.03

# density of only species 1
density(mybactgrid, 1)  # approx 0.01

# Making a contourplot
# --------------------

# simulate your data, in your case output from model
simulate(nbact, nphage) = log(nbact^2 + 2nphage + nbact*nphage - 0.1nbact^2) + 0.1randn()

# these will be a list of your number of species?
nbactvals = 1:10:1000
nphagevals = 1:100:10000

results = zeros(length(nphagevals), length(nbactvals))

# perform simulations

for (i, b) in enumerate(nphagevals)
    for (j, p) in enumerate(nbactvals)
        results[i, j] = simulate(b, p)
    end
end

using Plots

contour(nbactvals, nphagevals, results)
#contourf(nbactvals, nphagevals, results)  # use this to fill the contour
xlabel!("Number of bacteria")
ylabel!("Number of phages")
title!("My simulation")

# use these lines to plot in logaritmic scale
#xaxis!(:log)
#yaxis!(:log)
