#=
Created on 02/09/2021 14:37:03
Last update: Friday 03 September 2021

@author: Michiel Stock
michielfmstock@gmail.com

Analysis on the generalists dataset.
=#

using DrWatson
quickactivate(@__DIR__, "PhageSim")
using Statistics
using Plots
using LinearAlgebra
using CSV, DataFrames
using PhageSim

burnin = 400  # time to let the system settle
nspecies = 3
legend=true

summaries = DataFrame()

for θ in [0, 0.1, 0.2, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45,
                    0.475, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    results = CSV.read(datadir("generalists/generalists_$(θ).csv"), DataFrame)
    results[!,:entropy_bacteria] = results[:, ["count_bacteria_$i" for i in 1:2]] |> eachrow .|> Tuple .|> entropy
    results[!, :generalist] .= results[:, :count_phages_3] .> 0
    results[!, :specialists] .= (results[:, :count_phages_1] .> 0) .| (results[:, :count_phages_2] .> 0)
    results[!, :fgeneralists] .= results[:, :count_phages_3] ./ results[:, :count_phages]

    summary = results[results.step.≥burnin,:] |>
            (df -> groupby(df, :ensemble)) |>
            (df -> combine(df, :count_bacteria=>mean,
                                :count_phages=>mean,
                                :entropy_bacteria=>mean,
                                :generalist=>last,
                                :specialists=>last,
                                :fgeneralists=>mean))

    summary[:,:θ] .= θ

    append!(summaries, summary)
             

    reps = maximum(results.ensemble)
        for rep in 1:reps

            # plot growth curves
        
            results_rep = results[results.ensemble.==rep,:]

            pbact = plot(;legend)
            for i in 1:2
                plot!(pbact, results_rep[!, "count_bacteria_$i"], label="bacteria sp. $i", color=agentcolor(i))
            end
            xlabel!(pbact, "time")
            title!(pbact, "Bacteria ($nspecies species)")

            pphage = plot(;legend)
            for i in 1:3
                plot!(pphage, results_rep[!, "count_phages_$i"], label="phages sp. $i", color=agentcolor(i))
            end
            xlabel!(pphage, "time")
            title!(pphage, "Phages ($nspecies species)")

            plot(pbact, pphage, layout=(2, 1))

            savefig(plotsdir("generalists/growth_curves/growth_θ=$(θ)_rep=$rep.pdf"))

        end
end

summary = combine(groupby(summaries, :θ),
                :count_bacteria_mean => mean => :av_bact,
                :count_phages_mean => mean => :av_phages,
                :entropy_bacteria_mean => mean => :av_bact_entropy,
                :generalist_last => mean => :Pgen,
                :specialists_last => mean => :Pspec,
                :fgeneralists_mean  => mean => :fgeneralists)


CSV.write(datadir("generalists/summary.csv"), summary)

# make plot of the results

using Plots

p_fractions = plot(xlabel="θ", legend=:right, title="Average bacterial diversity", palette=:lightrainbow)
plot!(p_fractions, summary.θ, summary.Pspec, label="fraction of specialist", marker=:o, lw=2)
plot!(p_fractions, summary.θ, summary.Pgen, label="fraction of generalists", marker=:square, lw=2)
plot!(p_fractions, summary.θ, summary.av_bact_entropy, label="Shannon index", marker=:^, lw=2)

p_numbers = plot(xlabel="θ", legend=:left, title="Average number of agents",  palette=:lightrainbow)
plot!(p_numbers, summary.θ, summary.av_bact, label="average number of bacteria", lw=2, marker=:o)
plot!(p_numbers, summary.θ, summary.fgeneralists .* summary.av_bact, label="average number of bacterial generalists", lw=2, marker=:square)
plot!(p_numbers, summary.θ, summary.av_phages, label="average number of phages", lw=2, marker=:^)

plot(p_fractions, p_numbers, layout=(2,1), size=(600, 500))

savefig(plotsdir("generalists/summary_generalists.pdf"))