#=
Created on Friday 16 October 2020
Last update: Wednesday 31 March 2020

@author: Michiel Stock
michielfmstock@gmail.com

Assessing the effect of the structure of the interaction matrix on the bacteria-phage matrix.
File to plot the results.
=#

using DrWatson
quickactivate(@__DIR__, "PhageSim")
using PhageSim
using Statistics
using Plots, StatsPlots
using Agents, AgentsPlots
using LinearAlgebra
using CSV, DataFrames

for infecttype in ["uniform", "unique", "nested", "secundair"]
    results = CSV.read(datadir("interactions/$infecttype.csv"), DataFrame)

    #@df results corrplot([:count_bacteria1, :count_bacteria2, :count_bacteria3, :count_phages1, :count_phages2, :count_phages3])
    corrplot(results[:,5:end-1] |> Matrix, labels=names(results)[5:end-1], size=(1500,1500))
    savefig(plotsdir("interactions/$(infecttype)_corplot.png"))

    # plot growth curves
    for rep in 1:50
       
        results_rep = results[results.replicate.==rep,:]

        pbact = plot(results_rep.count_bacteria1, label="bacteria sp. 1", color=:red)
        plot!(results_rep.count_bacteria2, label="bacteria sp. 2", color=:blue)
        plot!(results_rep.count_bacteria3, label="bacteria sp. 3", color=:green)
        xlabel!("time")
        title!("Bacteria")

        pphage = plot(results_rep.count_phages1, label="phages sp. 1", color=:red, ls=:dash)
        plot!(results_rep.count_phages2, label="phages sp. 2", color=:blue, ls=:dash)
        plot!(results_rep.count_phages3, label="phages sp. 3", color=:green, ls=:dash)
        xlabel!("time")
        title!("Phages")

        plot(pbact, pphage, layout=(2, 1))

        savefig(plotsdir("interactions/growth_curves/growth_$(infecttype)_$rep.pdf"))
    end
end

