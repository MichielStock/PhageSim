#=
Created on 30/08/2021 17:12:56
Last update: 1 September 2021

@author: Michiel Stock
michielfmstock@gmail.com

Simulates the recovery of a system when there is only a single dominat species
=#


using DrWatson, Distributed
quickactivate(@__DIR__, "PhageSim")

using InteractiveDynamics, CairoMakie

addprocs(8)

@everywhere begin
    #using Pkg; Pkg.activate(".")
    using PhageSim
    using Agents
    using BSON, CSV, Statistics

    repl = 100
    tsteps = 500

    extent = (50, 50)

    nspecies = nbactsp = nphagesp = 10 

    nbact_sp1 = 400

    # structures of matrix
    p = 0.25  # set this depending on the number of sp
    if nspecies == 5
        p = 0.5
    elseif nspecies == 10
        p = 1.0
    end

    Pnone = 0.0
    Punif = p#/nspecies
    Punique = probs_unique(nspecies, p)
    Pnested = probs_nested(nspecies, p)
    Psec = probs_sec(nspecies, p, p/5)

    # general parameters
    burstsize = nspecies==10 ? 20.0 : 10.0 # set this depending on the number of sp, 20 for 10 sp, 10 for 3 sp
    ΔE = .2
    l = 0.5
    Δbact = l
    Δphage = Δbact
    pdie = 0.01
    pdecay = 0.1


    nbacteria = 100
    nphages = 0


    adata = [(bacteria, count), (phages, count)]

    # metaprogramming hocus pocus to get the counters
    for i in 1:nspecies
    push!(adata, (eval(Meta.parse("bacteria_$i")), count))
    push!(adata, (eval(Meta.parse("phages_$i")), count))
    end

    model_step_phages! = model_step_rand_phages(10)
end

         
for (infecttype, Pinf, Psym) in zip(["reference", "uniform", "unique", "nested", "secundair"],
                        [Pnone, Punif, Punique, Pnested, Psec],
                        [:Pnone, :Punif, :Punique, :Pnested, :Psec])

    println("Simulating $infecttype...")

    infecttype=="nested" || continue

    # hacky metaprogramming solution
    "@everywhere Pinf = $Psym" |> Meta.parse |> eval

    @everywhere function generator(seed)
        model = init_model(extent, min(extent...)/20; nbacteria, nphages, nbactsp, nphagesp,
                        burstsize, ΔE, l, Δbact, Δphage, pdie, pdecay, seed,
                        infection=infmodel(Pinf))
        # add additional bacteria of species 1
        maxiter = 100_000
        placed = 0
        for _ in 1:maxiter
            pos = (rand() * extent[1], rand() * extent[2])
            free = true
            for agent in nearby_agents(pos, model, 2l, exact=true)
                if isbacterium(agent)
                    free = false
                    break
                end
            end
            !free && continue
            # add agent of species 1
            add_agent!(pos, model, true, 1, rand())
            placed += 1
            placed ≥ nbact_sp1 && break
        end
        return model
    end

    parameters = @dict extent nbacteria nphages nbactsp nphagesp burstsize ΔE l Δbact Δphage pdie pdecay Pinf infecttype
    safesave(datadir("recovery/params_$(infecttype)_$nspecies.bson"), parameters)

    # no inititial phages, but add 10 randomly every step
    results, _, models = ensemblerun!(generator, agent_step!, model_step_phages!, tsteps; adata,
                                        ensemble=repl, parallel=true)

    anybact = [model.nbacteria > 0 for model in models] |> mean
    anyphage = [nagents(model) - model.nbacteria > 0 for model in models] |> mean

    println("Bacteria remain in $anybact of the models and phages remain in $anyphage of the models")

    safesave(datadir("recovery/$(infecttype)_$nspecies.csv"), results)

    # make simulation

    println("Making a simulation...")

    abm_video(
        plotsdir("recovery/animations/movie_$(infecttype)_$nspecies.mp4"),
        generator(1),
        agent_step!,
        model_step_phages!;
        ac=agentcolor,
        as=agentsize,
        title = "Model $infecttype",
        frames = tsteps,
        spf = 2,
        framerate = 5,
)
    println()
end