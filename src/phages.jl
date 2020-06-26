#=
Created on Friday 27 December 2019
Last update: Friday 26 June 2020

@author: Michiel Stock
michielfmstock@gmail.com

Functions to model the spread of the phages.
=#


export AbstractPhageRules, PhageRules, phagedecay, step_phages!

# PHAGE TYPE AND FUNCTIONS
# ------------------------

abstract type AbstractPhage <: AbstractAgent end

mutable struct Phage <: AbstractBacterium
    id::Int
    pos::Tuple{Int,Int}
    species::Int  # decribes the species of the phage
end

# RULES AND BEHAVIOUR
# -------------------


abstract type AbstractPhageRules end

struct PhageRules <: AbstractPhageRules
    pdecay::Float64
    R::Int
    function PhageRules(pdecay, R=1)
        @assert 0.0 ≤ pdecay ≤ 1.0 "`pdecay` should be in [0, 1], got $pdecay"
        @assert R > 0 "R should be a positive integer"
        return new(pdecay, R)
    end
end

"""
    decays(phage::AbstractPhage, phagerules::PhageRules)

Decay the phage?
"""
function decays(phage::AbstractPhage, phagerules::PhageRules)
    phagerules.pdecay == 0.0 && return false
    phagerules.pdecay == 1.0 && return true
    return phagerules.pdecay <= rand()
end


function agent_step!(phage::AbstractPhage, model)
    phagerules = pr(model.properties)  # get the phage rules
    interactionrules = ir(model.properties)  # get the interaction rules
    move!(phage, model)  # phages always move
    agents = get_node_agents(phage.pos, model)
    hosts = filter!(bacteria, agents)
    for bact in hosts
        # check if an infection occurs, immediately exit function if so
        infects!(phage, host, phagerules) && return
    end
    # decay phage
    decays(phage, phagerules) && kill_agent!(phage, model)
end

function move!(phage::AbstractPhage, model)
    phagerules = pr(model.properties)
    neighbors = node_neighbors(phage, model, r=phagerules.R)
    node = rand(neighbors)
    move_agent!(phage, node, model)
end
