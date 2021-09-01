#=
Created on 30/08/2021 16:08:09
Last update: Wednesday 1 September 2021

@author: Michiel Stock
michielfmstock@gmail.com

Experiment with the trade-off generalists vs specialists.
The system has two bacterial species and three phage species.
Each bacterium has one primary virus. The third phage is a generalist,

=#

using DrWatson, Distributed
quickactivate(@__DIR__, "PhageSim")

addprocs(8)

using InteractiveDynamics, CairoMakie

@everywhere begin
    using Pkg; Pkg.activate(".")
    using PhageSim, Agents
    using BSON, CSV, Statistics

    repl = 100
    tsteps = 500

    extent = (50, 50)

    # structures of matrix
    p = 0.25


    # general parameters
    burstsize = 10.0
    ΔE = .2
    l = 0.5
    Δbact = l
    Δphage = Δbact
    pdie = 0.01
    pdecay = 0.1


    nbacteria = 500
    nphages = 2000

    nbactsp = 2
    nphagesp = 3


    adata = [(bacteria, count), (phages, count),
            (bacteria_1, count), (bacteria_2, count),
            (phages_1, count), (phages_2, count), (phages_3, count)]
end

         
for θ in [0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1]

    println("Simulating θ=$θ...")

    "@everywhere θ = $θ" |> Meta.parse |> eval

    @everywhere Pinf = [p 0 θ*p;
            0 p θ*p;]

    @everywhere generator(seed) = init_model(extent, min(extent...)/20; nbacteria, nphages, nbactsp, nphagesp,
                        burstsize, ΔE, l, Δbact, Δphage, pdie, pdecay, seed,
                        infection=infmodel(Pinf))

    parameters = @dict extent nbacteria nphages nbactsp nphagesp burstsize ΔE l Δbact Δphage pdie pdecay Pinf
    safesave(datadir("generalists/params_generalists_$(θ).bson"), parameters)

    results, _, models = ensemblerun!(generator, agent_step!, model_step!, tsteps; adata,
                                    ensemble=repl, parallel=true)

    anybact = [model.nbacteria > 0 for model in models] |> mean
    anyphage = [nagents(model) - model.nbacteria > 0 for model in models] |> mean

    println("Bacteria remain in $anybact of the models and phages remain in $anyphage of the models")

    safesave(datadir("generalists/generalists_$(θ).csv"), results)

    # make simulation

    println("Making a simulation...")

    abm_video(
        plotsdir("generalists/movie_generalists_$θ.mp4"),
        generator(1),
        agent_step!,
        model_step!;
        ac=agentcolor,
        as=agentsize,
        title = "Model generalists θ=$θ",
        frames = tsteps,
        spf = 2,
        framerate = 5,
)

    println()
end