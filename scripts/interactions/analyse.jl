#=
Created on Friday 16 October 2020
Last update: Tuesday 5 October 2021

@author: Michiel Stock
michielfmstock@gmail.com

Assessing the effect of the structure of the interaction matrix on the bacteria-phage matrix.
File to plot the results.
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
add_total = true

summaries = DataFrame[]

for nspecies in [3, 5, 10]

    pphase = scatter(title="Phase plot ($nspecies species)",
                        xlabel="number of bacterial cells",
                        ylabel="number of phage particles", palette=:lightrainbow)

    #scatter!(pphase, [500], [1000], label="starting point")
    println("Loading and processing the results with $nspecies species...")

    for infecttype in ["reference", "uniform", "unique", "secundair", "nested"]
        results = CSV.read(datadir("interactions/$(infecttype)_$(nspecies).csv"), DataFrame)

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

        sample_to_plot = rand(1:reps)
        

        pphase_entropy = plot(title="Phase plot $infecttype ($nspecies species)",
                        xlabel="number of bacterial cells",
                        ylabel="diversity bacteria")
        for rep in 1:reps

            # plot growth curves
        
            results_rep = results[results.ensemble.==rep,:]

            pbact = plot(;legend)
            for i in 1:nspecies
                plot!(pbact, results_rep[!, "count_bacteria_$i"], label="bacteria sp. $i", color=agentcolor(i))
            end
            if add_total
                plot!(pbact, results_rep[!, "count_bacteria"], label="bacteria (total)", color="black", ls=:dash)
            end
            xlabel!(pbact, "time")
            title!(pbact, "Bacteria ($nspecies species)")

            pphage = plot(;legend)
            for i in 1:nspecies
                plot!(pphage, results_rep[!, "count_phages_$i"], label="phages sp. $i", color=agentcolor(i))
            end
            if add_total
                plot!(pphage, results_rep[!, "count_phages"], label="phages (total)", color="black", ls=:dash)
            end
            xlabel!(pphage, "time")
            title!(pphage, "Phages ($nspecies species)")

            plot(pbact, pphage, layout=(2, 1))

            savefig(plotsdir("interactions/growth_curves/$nspecies/growth_$(infecttype)_nsp=$(nspecies)_rep=$rep.pdf"))

            # plot diversity

            plot(xlabel="time", titel="Diversity ($nspecies species)", legend=:bottomleft)
            plot!(results_rep[!, :entropy_bacteria], label="entropy bacteria", lw=2)
            plot!(results_rep[!, :entropy_phages], label="entropy phages", lw=2)

            savefig(plotsdir("interactions/entropy/$nspecies/entropy_$(infecttype)_nsp=$(nspecies)_rep=$rep.pdf"))

            # make phase plots

            rep==sample_to_plot && plot!(pphase, results_rep[!, "count_bacteria"], results_rep[!, "count_phages"], label=infecttype, lw=2, alpha=0.8,)#, color="blue")
            #rep%10==0 && plot!(pphase_entropy, results_rep[!, "count_bacteria"], (2).^results_rep[!, :entropy_bacteria], alpha=0.9, label="")#, color="blue")
        end

        
        #savefig(pphase_entropy, plotsdir("interactions/phase_plots/$nspecies/phase_entropy_$(infecttype)_nsp=$(nspecies).pdf"))

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
    
    #scatter!(pphase, [500], [1000], label="starting point")
    savefig(pphase, plotsdir("interactions/phase_plots/phase_plot_nsp=$(nspecies).pdf"))

end


CSV.write(datadir("interactions/summaries.csv"), vcat(summaries...));
