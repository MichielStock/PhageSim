#=
Created on Monday 30 March 2020
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

This is a file for the bachelor student of the project
"DE IMPACT VAN BACTERIOFAGEN OP MICROBIËLE DIVERSITEIT".

Adapt the parameters and run this file. Your figures and simulation results
should be saved in the corresponding folders.
=#

using DrWatson
quickactivate(@__DIR__, "PhageSim")
using PhageSim
using Plots
import BSON: @save, @load

name = "mysim"  # give a name to your simulation

# PARAMETERS OF THE SIMULATION
# ----------------------------

# bacteria
# --------

prepr = 0.3
pmove = 0.4
pdie = 0.1

# phages
# -----

pdecay = 0.1  # probability that a viron decays
Rphages = 5  # mixing radius for the dispersion of the phages

# interactions
# ------------

# matrix to determine the chance of infection host-phage, nbact x nphages
Pinf =  [0.2 0.01 0.01; 0.01 0.2 0.01; 0.01 0.01 0.2]
burstsize = 10  # average number of virons
plysogeny = 0.0  # probability of entering the lysogentic cycle, either probability
# for all phages or a list of probabilities, one for each phage
plysis = 0.4  # probability that an infected bacterium lyses in every step
Rqs = 10  # radius for quorum sensing, lysis probability is proportional to the local density
# 0 means that it is independent

# simulation stuff
# ----------------

nsteps = 350
D = 100  # size of the grid
ninitbact = 1000
ninitphages = 4000

# PROCESSING THE PARAMETERS
# -------------------------

nbactsp, nphagesp = size(Pinf)

simpars = @dict nbactsp nphagesp prepr pmove pdie pdecay Rphages burstsize plysogeny plysis Rqs nsteps D ninitbact ninitphages

bactrules = BacteriaRules(prepr, pmove, pdie)
phagerules = PhageRules(pdecay, Rphages)
interactrules = InteractionRules(Pinf, burstsize, plysogeny, plysis, Rqs)

# SETTING UP THE INITIAL CONDITIONS
# ---------------------------------

# init the bacteria
bactgrid = BactGrid(nothing, D, D)  # empty grid
bactgrid[rand(CartesianIndices(bactgrid), ninitbact)] .= [Bacterium(rand(1:nbactsp)) for i in 1:ninitbact]

# init the phages
phagegrids = Array{Int,2}[]
for i in 1:nphagesp
    phagegrid = zeros(Int, D, D)
    phagegrid[rand(CartesianIndices(phagegrid), ninitphages)] .= 10  # place 10 phages
    push!(phagegrids, phagegrid)
end

# RUNNING THE SIMULATION
# ----------------------

# perform the simulations, the function `species_counts`, with keyword arguments
# `nspecies` for the number of bacterium species, count the number of bacteria
# and the number of phage sp in every time step

results = simulate!(bactgrid, phagegrids, bactrules, phagerules, interactrules,
                nsteps; resultfun=species_counts, nspecies=nbactsp)
