#=
Created on Saturday 28 December 2019
Last update: Thursday 02 July 2020

@author: Michiel Stock
michielfmstock@gmail.com

Implementation of the bacteria using Agents.jl
=#

export AbstractBacterium, Bacterium
export bacteria, prophage, haslatent, species, prophage!, density
export AbstractBacteriaRules, BacteriaRules, ProphageBacteriaRules,
            HeteroBacteriaRules, agent_step!

# BACTERIA TYPE AND FUNCTIONS
# ---------------------------

abstract type AbstractBacterium <: AbstractAgent end

mutable struct Bacterium <: AbstractBacterium
    id::Int
    pos::Tuple{Int,Int}
    species::Int  # decribes the species of the bacterium
    prophage::Union{Int,Nothing}  # either carries a latent phage (i) or not (0)
end

Bacterium(id, pos, species=1) = Bacterium(id, pos, species, nothing)

"""Test if agent is a bacterium."""
bacteria(a::AbstractAgent) = a isa AbstractBacterium

"""
    prophage(bact::AbstractBacterium)

Return the prophage of a bacterium. Return `nothing` if bact has no prophage.
"""
prophage(bact::AbstractBacterium) = bact.prophage

"""
    prophage!(bact::AbstractBacterium, sp=nothing)

Sets the prophage of `bact` to `sp`. Removes the prophage by default.
"""
prophage!(bact::AbstractBacterium, sp=nothing) = (bact.prophage = sp)


"""
    haslatent(bact::AbstractBacterium)

Test whether a bacterium has a laten prophage.
"""
haslatent(bact::AbstractBacterium) = !(bact.prophage isa Nothing)

haslatent(a::AbstractAgent) = false

"""
    species(bact::AbstractAgent)

Returns the species of a bacterium or phage.
"""
species(bact::AbstractAgent) = bact.species

"""
    species(bact::AbstractBacterium, sp::Int)

Test if `bact` is of species `sp`.
"""
species(bact::AbstractBacterium, sp::Int) = species(bact) == sp

"""
    density(bact, model, R, forspecies::Bool=true)

Compute the density of the bacteria in a region of radius `R` around `bact`.
If `species=false`, all bacteria will be taken into account. The default behaviour only
takes bacteria of the same species as `bact`.
"""
function density(bact, model, R, forspecies::Bool=true)
    #FIXME : does not work properly
    bactinregion = space_neighbors(bact.pos, model, R)
    nbact = forspecies ? count(b->species(b, species(bact)), bactinregion) : length(bactinregion)
    return nbact / (R + 1)^2
end


# RULES AND BEHAVIOUR
# -------------------

abstract type AbstractBacteriaRules end

"""
Rules for when all bacteria behave the same.
"""
struct BacteriaRules <: AbstractBacteriaRules
    prepr::Float64
    pmove::Float64
    pdie::Float64
    function BacteriaRules(prepr::Float64, pmove::Float64, pdie::Float64)
        @assert +(prepr, pmove, pdie) ≤ 1.0 &&
                all((prepr, pmove, pdie) .≥ 0) "behaviour of the bacteria should be valid probabilites"
        new(prepr, pmove, pdie)
    end
end

"""
Rules for when species of bacteria might show different behaviours.
"""
struct HeteroBacteriaRules <: AbstractBacteriaRules
    prepr::Array{Float64,1}
    pmove::Array{Float64,1}
    pdie::Array{Float64,1}
end

function BacteriaRules(prepr::Array{Float64,1}, pmove::Array{Float64,1}, pdie::Array{Float64,1})
    @assert all(.+(prepr, pmove, pdie) .≤ 1.0) && all(prepr .≥ 0) && all(pmove .≥ 0) &&
             all(pdie .≥ 0) "behaviour of the bacteria should be valid probabilites"
    HeteroBacteriaRules(prepr, pmove, pdie)
end

"""
Rules for when a bacterium with a prophage exhibits a different behaviour from
a bacterium without a prophage.
"""
struct ProphageBacteriaRules <: AbstractBacteriaRules
    probsnoinf::Tuple{Float64,Float64,Float64}
    probsinf::Tuple{Float64,Float64,Float64}
    precover::Float64
end

function BacteriaRules(probsnoinf::Tuple{Float64,Float64,Float64},
                probsinf::Tuple{Float64,Float64,Float64},
                precover::Float64=0.0)
    @assert sum(probsnoinf) ≤ 1.0 &&
            all(probsnoinf .≥ 0) "behaviour of the bacteria should be valid probabilites"
    @assert sum(probsinf) ≤ 1.0 &&
            all(precover .≥ 0) "behaviour of the bacteria should be valid probabilites"
    @assert 0 ≤ precover ≤ 1 "`precover` should be a valid probability"
    ProphageBacteriaRules(probsnoinf, probsinf, precover)
end

"""
    bacteriaprobs(br::BacteriaRules, bact::AbstractBacterium)

Get the behaviour parameters for a specific bacterium.
"""
bacteriaprobs(bact::AbstractBacterium, br::BacteriaRules) = (br.prepr, br.pmove, br.pdie)

"""
    bacteriaprobs(br::ProphageBacteriaRules, bact::AbstractBacterium)

Get the behaviour parameters depending on whether the bacterium has a prophage
or not.
"""
bacteriaprobs(bact::AbstractBacterium, br::HeteroBacteriaRules) =
                    species(bact) |> i -> (br.prepr[i], br.pmove[i], br.pdie[i])

bacteriaprobs(bact::AbstractBacterium, br::ProphageBacteriaRules) =
        haslatent(bact) ? br.probsinf : br.probsinf

"""
By default an infected phage cannot recover.
"""
recovers(bact::AbstractBacterium, br::AbstractBacteriaRules) = false

"""
Computes the probability that a bacterium with a prophage looses recovers.
"""
recovers(bact::AbstractBacterium, br::ProphageBacteriaRules) = br.precover > 0.0 && rand() < br.precover


"""
    agent_step!(bact::Bacterium, model)

Defines a step of a bacterium.
"""
function agent_step!(bact::AbstractBacterium, model)
    bactrules = br(model.properties)  # get the bacteria rules
    # if the bacterium has a prophage it might lyse or recover
    if haslatent(bact)
        if recovers(bact, bactrules)
            bact.prophage = nothing  # cured
        else
            interactionrules = ir(model.properties)
            lyses(bact, interactionrules, model) && return lyse!(bact, model)
        end
    end
    # get behaviour parameters
    prepr, pmove, pdie = bacteriaprobs(bact, bactrules)
    r = rand()
    if r < pmove
        return move!(bact, model)
    elseif r < pmove + prepr
        return reproduce!(bact, model)
    elseif r < pmove + prepr + pdie
        return kill_agent!(bact, model)
    end
end

"""Check if there are no bacteria on a node."""
isfree(node, model) = !any(bacteria, get_node_agents(node, model))

function move!(bact::AbstractBacterium, model)
    neighbors = node_neighbors(bact, model)
    node = rand(neighbors)
    if isfree(node, model)
        move_agent!(bact, node, model)
    end
end

function reproduce!(bact::AbstractBacterium, model)
    neighbors = node_neighbors(bact, model)
    node = rand(neighbors)
    if isfree(node, model)
        id = nextid(model)
        add_agent!(Bacterium(id, node, bact.species, bact.prophage), model)
    end
end
