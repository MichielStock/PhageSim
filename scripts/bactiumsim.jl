#=
Created on Friday 31 January 2020
Last update: Friday 10 July 2020

@author: Michiel Stock
michielfmstock@gmail.com

Simple simulation of the bacteria without the phages.
=#

using DrWatson
quickactivate(@__DIR__, "PhageSim")
using PhageSim, Plots
using Agents, AgentsPlots

tsteps = 100
N = 100
ninitbact = 25
nspecies = 1
pmove = 0.2
Ediv = 2.0
Emax = 10.0
ΔE = 0.2
σE = 0.1

space = GridSpace((N, N), moore=true)

bactrules = BacteriaRules(pmove, Ediv, ConstantEnergyUpdate(ΔE, Emax))

model = AgentBasedModel(Bacterium, space, scheduler=random_activation,
                        properties=bactrules)

for _ in 1:ninitbact
    id = nextid(model)
    species = rand(1:nspecies)
    pos = (rand(1:N), rand(1:N))
    bact = Bacterium(id, pos, species)
    add_agent!(bact, model)
    move_agent_single!(bact, model)
end

adata = [(bacteria, count)]

results, _ = run!(model, agent_step!, tsteps, adata=adata)

plot(results.count_bacteria)
