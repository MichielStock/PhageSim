#alle figuren voor onderzoeksvraag1
#worden in 1 keer aangemaakt

using DrWatson
quickactivate(@__DIR__, "PhageSim")
using PhageSim
using Plots
import BSON: @save, @load

#scenario
#ninitbact, pdie, prepr
ldens = [50, 0.15, 0.30]
mdens = [200, 0.10, 0.35]
hdens = [500, 0.05, 0.40]
dens = [ldens, mdens, hdens]

#situatie
#plysogeny, plysis
sit1 = [0.0, 0.5]
sit2 = [0.0, 0.1]
sit3 = [0.1, 0.3]
sit4 = [0.3, 0.6]
sit5 = [0.5, 0.5]
sit6 = [0.5, 0.3]
sit = [sit1, sit2, sit3, sit4, sit5, sit6]

for i = 1:length(dens)
    for j = 1:length(sit)

        name = "scenario" * string(i) * "situatie" * string(j)  # give a name to your simulation

        # PARAMETERS OF THE SIMULATION
        # ----------------------------

        # bacteria
        # --------

        prepr = dens[i][3]
        pmove = 0.4
        pdie = dens[i][2]

        # phages
        # -----

        pdecay = 0.01  # probability that a viron decays
        Rphages = 5  # mixing radius for the dispersion of the phages

        # interactions
        # ------------

        # matrix to determine the chance of infection host-phage, nbact x nphages
        Pinf =  [0.2 0.01 0.01; 0.01 0.2 0.01; 0.01 0.01 0.2]
        burstsize = 10  # average number of virons
        plysogeny = sit[j][1]  # probability of entering the lysogentic cycle, either probability
        # for all phages or a list of probabilities, one for each phage
        plysis = sit[j][2]  # probability that an infected bacterium lyses in every step
        Rqs = 0  # radius for quorum sensing, lysis probability is proportional to the local density
        # 0 means that it is independent

        # simulation stuff
        # ----------------

        nsteps = 200
        D = 100  # size of the grid
        ninitbact = floor(Int64, dens[i][1])
        ninitphages = 100

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

        # DATA STORAGE AND VISUALIZATION
        # ------------------------------

        # generate a unique name
        fname = savename(name, simpars)

        # save the simulation in a BSON file
        @save datadir(fname) * ".bson" results

        # can be loaded using
        # @load datadir(fname) * ".bson"

        # plotting

        bactres = [res[1][sp] for res in results, sp in 1:nbactsp]
        phageres = [res[2][sp] for res in results, sp in 1:nphagesp]

        p = plot(plot(bactres, labels=["sp. 1" "sp. 2" "sp. 3"], xlabel="step",
                ylabel="number of bacteria", title="Bacteria composition"),
            plot(phageres, labels=["sp. 1" "sp. 2" "sp. 3"], ls=:dash, xlabel="step",
                    ylabel="number of phages", title="Phage composition"))

        savefig(p, plotsdir(fname) * ".png")
    end
end

a = 3
