#=
Created on Saturday 28 December 2019
Last update: Thursday 02 July 2020

@author: Michiel Stock
michielfmstock@gmail.com

Implementation of the bacteria using Agents.jl
=#

export AbstractBacterium, Bacterium
export bacteria, prophage, haslatent, species, prophage!, density, energy, energy!
export AbstractBacteriaRules, BacteriaRules
export AbstractEnergyupdate, ConstantEnergyUpdate, RandomEnergyUpdate

# BACTERIA TYPE AND FUNCTIONS
# ---------------------------

abstract type AbstractBacterium <: AbstractAgent end

mutable struct Bacterium <: AbstractBacterium
    id::Int
    pos::Tuple{Int,Int}
    species::Int  # decribes the species of the bacterium
    energy::Float64
    prophage::Union{Int,Nothing}  # either carries a latent phage (i) or not (0)
end

Bacterium(id, pos, species=1, energy=1.0) = Bacterium(id, pos, species, energy, nothing)

"""Test if agent is a bacterium."""
bacteria(a::AbstractAgent) = a isa AbstractBacterium

"""
    energy(bact::AbstractBacterium)

Get the energy level of a bacterium.
"""
energy(bact::AbstractBacterium) = bact.energy

"""
    energy!(bact::AbstractBacterium, level::Float64)

Set the energy level of a bacterium.
"""
energy!(bact::AbstractBacterium, level::Float64) = (bact.energy = level)


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
    bactinregion = space_neighbors(bact.pos, model, R)
	nbact = count(a->bacteria(a) && (!forspecies || species(a)==species(bact)),
							model.agents[bactinregion])
    return nbact / (R + 1)^2
end


# RULES AND BEHAVIOUR
# -------------------

abstract type AbstractEnergyupdate end

struct ConstantEnergyUpdate <: AbstractEnergyupdate
	ΔE::Float64
	Emax::Float64
end

struct RandomEnergyUpdate <: AbstractEnergyupdate
	ΔE::Float64
	Emax::Float64
	σE::Float64
end

function energyupdate(bact::AbstractBacterium, ur::ConstantEnergyUpdate)
	bact.energy += ur.ΔE
	bact.energy > ur.Emax && energy!(bact, ur.Emax)
end

function energyupdate(bact::AbstractBacterium, ur::RandomEnergyUpdate)
	bact.energy += ur.ΔE + ur.σE * randn()
	bact.energy > ur.Emax && energy!(bact, ur.Emax)
end

abstract type AbstractBacteriaRules end

struct BacteriaRules{TPM,TED,TUF<:AbstractEnergyupdate} <: AbstractBacteriaRules
    pmove::TPM
	Ediv::TED
	Eupdate::TUF
end

br(bactrules::BacteriaRules) = bactrules

#TODO write constructor for this stuff

getspeciespar(bact::AbstractBacterium, par::Number) = par
getspeciespar(bact::AbstractBacterium, par::AbstractVector) = par[species(bact)]

pmove(bact, br) = getspeciespar(bact, br.pmove)
Emax(bact, br) = getspeciespar(bact, br.Emax)
Ediv(bact, br) = getspeciespar(bact, br.Ediv)
ΔEmain(bact, br) = getspeciespar(bact, br.ΔEmain)
Eupdate(br) = br.Eupdate

"""
    agent_step!(bact::Bacterium, model)

Defines a step of a bacterium.
"""
function agent_step!(bact::AbstractBacterium, model)
    bactrules = br(model.properties)  # get the bacteria rules
    # if the bacterium has a prophage it might lyse or recover
    if haslatent(bact)
        interactionrules = ir(model.properties)
        lyses(bact, interactionrules, model) && return lyse!(bact, model)
    end
	# update energy
	energyupdate(bact, Eupdate(bactrules))
	energy(bact) < 0.0 && return kill_agent!(bact, model)
	# try to move
    pmove(bact, bactrules) ≤ rand() && move!(bact, model)
    (energy(bact) > Ediv(bact, bactrules)) && reproduce!(bact, model)
end

"""Check if there are no bacteria on a node."""
isfree(pos, model) = !any(bacteria, get_node_agents(pos, model))

function move!(bact::AbstractBacterium, model)
    neighbors = node_neighbors(bact, model)
    pos = rand(neighbors)
    if isfree(pos, model)
        move_agent!(bact, pos, model)
    end
end

function reproduce!(bact::AbstractBacterium, model)
    neighbors = node_neighbors(bact, model)
    pos = rand(neighbors)
    if isfree(pos, model)
		E = energy(bact) / 2
		energy!(bact, E)
        id = nextid(model)
        add_agent!(Bacterium(id, pos, bact.species, E, bact.prophage), pos, model)
    end
end
