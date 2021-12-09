#=
Created on Thursday 2 September 2021
Last update: Thursday 2 September 2021

@author: Michiel Stock
michielfmstock@gmail.com

Analysing the effect of recovery
=#

using DrWatson
quickactivate(@__DIR__, "PhageSim")
using Statistics
using Plots, TernaryPlots
using LinearAlgebra
using CSV, DataFrames
using PhageSim: agentcolor, entropy

legend = false
burnin = 200  # time to let the system settle

summaries = DataFrame[]

for nspecies in [3, 5, 10]
    println("Loading and processing the results with $nspecies species...")

    for infecttype in ["reference", "uniform", "unique", "nested", "secundair"]
        results = CSV.read(datadir("recovery/$(infecttype)_$(nspecies).csv"), DataFrame)

        infecttype = infecttype=="secundair" ? "secondary" : infecttype
        infecttype = infecttype=="reference" ? "control" : infecttype

        results[!,:entropy_bacteria] = results[:, ["count_bacteria_$i" for i in 1:nspecies]] |> eachrow .|> Tuple .|> entropy
        results[!,:entropy_phages] = results[:, ["count_phages_$i" for i in 1:nspecies]] |> eachrow .|> Tuple .|> entropy

        results[!,:bacteria_sp] = results[:, ["count_bacteria_$i" for i in 1:nspecies]] |> eachrow .|> Tuple .|> (t->count(>(0), t))
        results[!,:phage_sp] = results[:, ["count_phages_$i" for i in 1:nspecies]] |> eachrow .|> Tuple .|> (t->count(>(0), t))

        #@df results corrplot([:count_bacteria1, :count_bacteria2, :count_bacteria3, :count_phages1, :count_phages2, :count_phages3])
        #corrplot(results[:,5:end-1] |> Matrix, labels=names(results)[5:end-1], size=(1500,1500))
        #savefig(plotsdir("interactions/$(infecttype)_$(nspecies)_corplot.png"))

        

        reps = maximum(results.ensemble)
        for rep in 1:reps

            # plot growth curves
        
            results_rep = results[results.ensemble.==rep,:]

            pbact = plot(;legend)
            for i in 1:nspecies
                plot!(pbact, results_rep[!, "count_bacteria_$i"], label="bacteria sp. $i", color=agentcolor(i))
            end
            xlabel!(pbact, "time")
            title!(pbact, "Bacteria ($nspecies species)")

            pphage = plot(;legend)
            for i in 1:nspecies
                plot!(pphage, results_rep[!, "count_phages_$i"], label="phages sp. $i", color=agentcolor(i))
            end
            xlabel!(pphage, "time")
            title!(pphage, "Phages ($nspecies species)")

            plot(pbact, pphage, layout=(2, 1))

            savefig(plotsdir("recovery/growth_curves/$nspecies/growth_$(infecttype)_nsp=$(nspecies)_rep=$rep.pdf"))

            # plot diversity

            plot(xlabel="time", titel="Diversity ($nspecies species)", legend=:bottomleft)
            plot!(results_rep[!, :entropy_bacteria], label="entropy bacteria", lw=2)
            plot!(results_rep[!, :entropy_phages], label="entropy phages", lw=2)

            savefig(plotsdir("recovery/entropy/$nspecies/entropy_$(infecttype)_nsp=$(nspecies)_rep=$rep.pdf"))
        end

        # check entropy, density etc

        

        summary = results |> (df -> df[df.step.≥burnin,:]) |>
                    (df -> groupby(df, :ensemble)) |>
                    (df -> combine(df, :count_bacteria => mean => :number_bacteria,
                                        :count_bacteria => std => :std_bacteria,
                                        :count_phages => mean =>:number_phages,
                                        :count_phages => std => :std_phages,
                                        [:count_bacteria, :count_phages] => cor => :cor_bact_phages,
                                        :bacteria_sp => last => :final_bact_sp,
                                        :phage_sp => last => :final_phage_sp,
                                        :entropy_bacteria => mean => :entropy_bacteria,
                                        :entropy_phages => mean => :entropy_phages))
        summary[!,:nspecies] .= nspecies
        summary[!, :infecttype] .= infecttype

        push!(summaries, summary)
    end

end


CSV.write(datadir("recovery/summaries.csv"), vcat(summaries...))
