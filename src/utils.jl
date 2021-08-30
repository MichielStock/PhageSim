#=
Created on Wednesday 26 March 2020
Last update: Thurday 26 August 2021

@author: Michiel Stock
michielfmstock@gmail.com

General utilies for modelling.
=#

using LinearAlgebra
export agentcolor, agentsize, bacteria, phages
export probs_nested, probs_sec, probs_unique

agentcolor(a::Agent) = agentcolor(species(a))
agentcolor(i::Integer) = [:red, :blue, :green, :yellow, :black,
                    :pink, :orange, :black, :purple, :grey][i]
                    
agentsize(agent::Agent) = isbacterium(agent) ? 15.0 : 2.0
agentmarker(agent) = isbacterium(agent) ? :circle : :^

bacteria(a) = isbacterium(a)
phages(a) = !isbacterium(a)

"""Create interaction matrix for `nspecies`, where each species has
one partner with probability `p` and all other partner with probability `psec`."""
probs_sec(nspecies, p, psec) = [i==j ? p : psec
                        for i in 1:nspecies, j in 1:nspecies]

"""Create interaction matrix for `nspecies`, where each species has
one partner with probability `p`."""
probs_unique(nspecies, p) = Matrix(p * I, nspecies, nspecies)

"""Create a nested interaction matrix for `nspecies`with base probability `p`."""
probs_nested(nspecies, p) = [(i ≥ j) * p / min(i, nspecies-j+1) for i in 1:nspecies, j in 1:nspecies]


for i in 1:20
    "bacteria_$i(a::Agent) = isbacterium(a) && species(a) == $i" |> Meta.parse |> eval
    "phages_$i(a::Agent) = !isbacterium(a) && species(a) == $i" |> Meta.parse |> eval
    "export phages_$i, bacteria_$i" |> Meta.parse |> eval
end


