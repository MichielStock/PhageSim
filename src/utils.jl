#=
Created on Wednesday 26 March 2020
Last update: Thurday 26 August 2021

@author: Michiel Stock
michielfmstock@gmail.com

General utilies for modelling.
=#

export agentcolor, agentsize

agentcolor(agent) = [:red, :blue, :green, :yellow, :black, :pink][agent.species]
agentsize(agent) = isbacterium(agent) ? 15.0 : 2.0
agentmarker(agent) = isbacterium(agent) ? :circle : :^


# TODO: export
for i in 1:20
    "bacteria_$i(a::Agent) = isbacterium(a) && species(a) == $i" |> Meta.parse |> eval
    "phage_$i(a::Agent) = !isbacterium(a) && species(a) == $i" |> Meta.parse |> eval
end


