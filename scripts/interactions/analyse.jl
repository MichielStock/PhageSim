#=
Created on Friday 16 October 2020
Last update: Monday 30 August 2021

@author: Michiel Stock
michielfmstock@gmail.com

Assessing the effect of the structure of the interaction matrix on the bacteria-phage matrix.
File to plot the results.
=#

using DrWatson
quickactivate(@__DIR__, "PhageSim")
using Statistics
using Plots#, StatsPlots
using LinearAlgebra
using CSV, DataFrames
using PhageSim: agentcolor

nspecies = 3
legend = false

for infecttype in ["reference", "uniform", "unique", "nested", "secundair"]
    results = CSV.read(datadir("interactions/$(infecttype)_$(nspecies).csv"), DataFrame)

    #@df results corrplot([:count_bacteria1, :count_bacteria2, :count_bacteria3, :count_phages1, :count_phages2, :count_phages3])
    #corrplot(results[:,5:end-1] |> Matrix, labels=names(results)[5:end-1], size=(1500,1500))
    #savefig(plotsdir("interactions/$(infecttype)_$(nspecies)_corplot.png"))

    # plot growth curves

    reps = maximum(results.ensemble)
    for rep in 1:reps
       
        results_rep = results[results.ensemble.==rep,:]

        pbact = plot(;legend)
        for i in 1:nspecies
            plot!(pbact, results_rep[!, "count_bacteria_$i"], label="bacteria sp. $i", color=agentcolor(i))
        end
        xlabel!(pbact, "time")
        title!(pbact, "Bacteria")

        pphage = plot(;legend)
        for i in 1:nspecies
            plot!(pphage, results_rep[!, "count_phages_$i"], label="phages sp. $i", color=agentcolor(i))
        end
        xlabel!(pphage, "time")
        title!(pphage, "Phages")

        plot(pbact, pphage, layout=(2, 1))

        savefig(plotsdir("interactions/growth_curves/growth_$(infecttype)_nsp=$(nspecies)_rep=$rep.pdf"))
    end
end

