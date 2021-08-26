#=
Created on 31/03/2021 09:58:43
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

When will lysogeny dominate?
=#

using DrWatson
quickactivate(@__DIR__, "PhageSim")
using PhageSim, Plots
using Agents, AgentsPlots
using LinearAlgebra
using BSON, CSV

repl = 50
tsteps = 1000

dims = (100, 100)

# structures of matrix
p = 0.1

# general parameters
burstsize = 5.0
pdecay = 0.1
plysogeny = [0.0, 0.25, 0.5, 0.75,1.0]
plysis = 0.05
pmovebact = 0.5
ΔE = 2.0
Ediv = 10.0
Emax = 25.0
σE = 2.0

n_bacteria =1000
n_phages = 25_000

parameters = Dict(
    :burstsize => 1.0,
    :pdecay => 0.1,
    :plysogeny => plysogeny,
    :plysis => 0.0,
    :pmovebact => 0.5,
    :ΔE => 2.0,
    :Ediv => 10.0,
)

adata = [(bacteria, count), (phages, count), (haslatent, count),
        (bacteria1, count), (bacteria2, count), (bacteria3, count), (bacteria4, count), (bacteria5, count),
        (phages1, count), (phages2, count), (phages3, count), (phages4, count), (phages5, count)]




for ΔE in [1.0, 2.0, 5.0, 10.0] 

    println("Simulating ΔE=$ΔE...")

    parameters[:ΔE] = ΔE

    safesave(datadir("lysogeny/parameters_ΔE=$ΔE.bson"), parameters)

    # generate basic rules
    rules = PhageSimRules(Pinf=p*ones(5,5), burstsize=burstsize, pdecay=pdecay, plysogeny=plysogeny,
                    plysis=plysis, pmovebact=pmovebact, Ediv=Ediv,
                    Eupdate=RandomEnergyUpdate(ΔE, Emax, σE))
    
    model = initialize_model(rules,
                    n_bacteria=n_bacteria,
                    n_phages=n_phages,
                    n_phage_sp=5,
                    n_bact_sp=5,
                    dims=dims)


    results, _ = run!(model, agent_step!, tsteps, adata=adata, replicates=repl, parallel=true)
    safesave(datadir("lysogeny/results_ΔE=$ΔE.csv"), results)

end