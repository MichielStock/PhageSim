#=
Created on 31/03/2021 10:52:59
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Analyse the results of exploring lysogeny
=#

using DrWatson
quickactivate(@__DIR__, "PhageSim")
using PhageSim
using Statistics
using Plots, StatsPlots
using Agents, AgentsPlots
using LinearAlgebra
using CSV, DataFrames

for ΔE in [1.0, 2.0, 5.0, 10.0] 
    results = CSV.read(datadir("lysogeny/results_ΔE=$ΔE.csv"), DataFrame)

    #@df results corrplot([:count_bacteria1, :count_bacteria2, :count_bacteria3, :count_phages1, :count_phages2, :count_phages3])
    #corrplot(results[:,5:end-1] |> Matrix, labels=names(results)[5:end-1], size=(2500,2500))
    #savefig(plotsdir("lysogeny/ΔE=$(ΔE)_corplot.png"))

    nrep = maximum(results.replicate)
    # plot growth curves
    for rep in 1:nrep
       
        results_rep = results[results.replicate.==rep,:]

        pbact = plot(results_rep.count_bacteria, label="bacteria", color=:black)
        xlabel!("time")
        title!("Bacteria")

        pphage = plot(results_rep.count_phages1, label="phages sp. 1 (0.0)", color=:red, ls=:dash)
        plot!(results_rep.count_phages2, label="phages sp. 2 (0.25)", color=:blue, ls=:dash)
        plot!(results_rep.count_phages3, label="phages sp. 3 (0.5)", color=:green, ls=:dash)
        plot!(results_rep.count_phages4, label="phages sp. 4 (0.75)", color=:yellow, ls=:dash)
        plot!(results_rep.count_phages5, label="phages sp. 5 (1.0)", color=:pink, ls=:dash)
        xlabel!("time")
        title!("Phages")

        plot(pbact, pphage, layout=(2, 1))

        savefig(plotsdir("lysogeny/growth_curves/growth_curve_ΔE=$(ΔE)_$rep.pdf"))
    end
end
