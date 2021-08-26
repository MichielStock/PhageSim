#=
Created on Thursday 16 July 2020
Last update: Friday 09 October 2020

@author: Michiel Stock
michielfmstock@gmail.com

Test of simple host-phage system.
=#

using DrWatson
quickactivate(@__DIR__, "PhageSim")
using PhageSim, Plots
using Agents, AgentsPlots
using BSON, CSV


sim_name = "specific"
repl=1
tsteps = 1000

dims = (50, 50)


Pinf = [0.3 0 0; 0 0.3 0; 0 0 0.3]

burstsize = 1.0
pdecay = 0.05
plysogeny = 0.0
plysis = 0.05
pmovebact = 0.5
ΔE = 2.0
Ediv = 10.0

# generate basic rules
rules = PhageSimRules(Pinf=Pinf, burstsize=burstsize, pdecay=pdecay, plysogeny=plysogeny,
                plysis=plysis, pmovebact=pmovebact, Ediv=Ediv,
                Eupdate=RandomEnergyUpdate(ΔE, 20.0, 1.0))

model = initialize_model(rules,
                        n_bacteria=500,
                        n_phages=10_000,
                        n_phage_sp=3,
                        n_bact_sp=3,
                        dims=dims)

adata = [(bacteria, count), (phages, count), (haslatent, count),
        (bacteria1, count), (bacteria2, count), (bacteria3, count),
        (phages1, count), (phages2, count), (phages3, count)]

results, _ = run!(model, agent_step!, tsteps, adata=adata, replicates=repl, parallel=true)

CSV.write(datadir(joinpath("results", "$sim_name.csv")), results)

results = results[results.replicate.==1,:]

pbact = plot(results.count_bacteria1, label="bacteria sp. 1", color=:red)
plot!(results.count_bacteria2, label="bacteria sp. 2", color=:blue)
plot!(results.count_bacteria3, label="bacteria sp. 3", color=:green)
xlabel!("time")
title!("Bacteria")

pphage = plot(results.count_phages1, label="phages sp. 1", color=:red, ls=:dash)
plot!(results.count_phages2, label="phages sp. 2", color=:blue, ls=:dash)
plot!(results.count_phages3, label="phages sp. 3", color=:green, ls=:dash)
xlabel!("time")
title!("Phages")

plot(pbact, pphage, layout=(2, 1))

savefig(plotsdir("$sim_name.pdf"))

parameters = @dict Pinf burstsize pdecay plysogeny plysis pmovebact ΔE Ediv
BSON.@save datadir(joinpath("parameters", "$sim_name.bson")) parameters results



#=
# making an animation
# -------------------

model = initialize_model(rules,
                        n_bacteria=300,
                        n_phages=300,
                        n_phage_sp=3,
                        n_bact_sp=3,
                        dims=(100, 100))

anim = @animate for t=0:500
        t > 0 && step!(model, agent_step!)
        plothp(model)
        title!("step $t")
end

gif(anim, plotsdir("threesp_test.gif"), fps = 5)
=#