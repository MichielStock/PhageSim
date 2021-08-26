#=
Created on Friday 09 October 2020
Last update: wedneday 31 March 2021

@author: Michiel Stock
michielfmstock@gmail.com

Assessing the effect of the structure of the interaction matrix on the bacteria-phage matrix.
=#

using DrWatson
quickactivate(@__DIR__, "PhageSim")
using PhageSim, Plots
using Agents, AgentsPlots
using LinearAlgebra
using BSON, CSV

repl = 5
tsteps = 1000

dims = (100, 100)

# structures of matrix
p = 0.1

Punif = p
Punique = 3p * Matrix(I, 3, 3)
Pnested = [p/3 0   0;
           p/3 p/2 0;
           p/3 p/2 p]
Psec = [p   0   p/2;
        p/2 p   0;
        0   p/2 p]

# general parameters
burstsize = 5.0
pdecay = 0.05
plysogeny = 0.0
plysis = 0.05
pmovebact = 0.5
ΔE = 2.0
Ediv = 10.0
Emax = 25.0
σE = 2.0

n_bacteria = 500
n_phages = 10_000

parameters = Dict{Symbol,Any}(
    :burstsize => 1.0,
    :pdecay => 0.1,
    :plysogeny => false,
    :plysis => 0.0,
    :pmovebact => 0.5,
    :ΔE => 2.0,
    :Ediv => 10.0,
)

adata = [(bacteria, count), (phages, count), (haslatent, count),
        (bacteria1, count), (bacteria2, count), (bacteria3, count),
        (phages1, count), (phages2, count), (phages3, count)]

         
for (infecttype, Pinf) in zip(["uniform", "unique", "nested", "secundair"],
                        [Punif, Punique, Pnested, Psec])

    println("Simulating $infecttype...")

    parameters[:Pinf] = Pinf

    safesave(datadir("interactions/$infecttype.bson"), parameters)

    # generate basic rules
    rules = PhageSimRules(Pinf=Pinf, burstsize=burstsize, pdecay=pdecay, plysogeny=plysogeny,
                    plysis=plysis, pmovebact=pmovebact, Ediv=Ediv,
                    Eupdate=RandomEnergyUpdate(ΔE, Emax, σE))
    
    model = initialize_model(rules,
                    n_bacteria=n_bacteria,
                    n_phages=n_phages,
                    n_phage_sp=3,
                    n_bact_sp=3,
                    dims=dims)


    results, _ = run!(model, agent_step!, tsteps, adata=adata, replicates=repl, parallel=true)
    safesave(datadir("interactions/$infecttype.csv"), results)

end