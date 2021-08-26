#=
Created on Friday 21 Feb 2020
Last update: Thursday 02 July 2020

@author: Michiel Stock
michielfmstock@gmail.com

Model the interactions between bacteria and their phages
=#


export AbstractInteractionRules, InteractionRules
export lysogenic, lyses, infects
using Distributions: Poisson

abstract type AbstractInteractionRules <: AbstractRules end

"""
Rules for the interactions between bacteria and phages.
"""
struct InteractionRules{TI,TB,TLG} <: AbstractInteractionRules
    Pinf::TI  # matrix to determine the chance of infection host-phage
    burstsize::TB  # average number of virons per unit of energy
    plysogeny::TLG  # probability of entering the lysogentic cycle
    plysis::Float64  # probability that an infected bacterium lyses
    Rqs::Int  # radius for quorum sensing
    function InteractionRules(Pinf, burstsize::Number,
            plysogeny::Union{BP,AbstractVector{BP}}=false,
            plysis=1.0, Rqs::Int=0)  where {BP <: Union{Bool, AbstractFloat}}
        @assert all(0 .≤ Pinf .≤ 1) "Pinf is Boolean or probabilities"
        @assert all(0 .≤ plysogeny .≤ 1) "All values for `plysogeny` should be Booleans or probabilities"
        @assert !(plysogeny isa AbstractVector) || size(Pinf, 2) == length(plysogeny) "`plysogeny` has incorrect size"
        @assert 0 ≤ plysis ≤ 1 "`plysis` should be a valid probability"
        @assert Rqs ≥ 0 "radius for quorum sensing `R` should be nonnegative"
        new{typeof(Pinf),typeof(burstsize),typeof(plysogeny)}(Pinf,burstsize,plysogeny,plysis,Rqs)
    end
end

"""Get the type of `Pinf`, used for dispatch."""
infectiontype(ir::InteractionRules{TI,TB,TLG}) where {TI<:Any,TB<:Any,TLG<:Any} = TI

"""Get the type of `burstsize`, used for dispatch."""
burssizetype(ir::InteractionRules{TI,TB,TLG}) where {TI<:Any,TB<:Any,TLG<:Any} = TB

"""Get the type of `plysogeny`, used for dispatch."""
plysogenytype(ir::InteractionRules{TI,TB,TLG}) where {TI<:Any,TB<:Any,TLG<:Any} = TLG

_lysogenic(phage, bact, interactionrules, ::Type{<:Bool}) = interactionrules.plysogeny
_lysogenic(phage, bact, interactionrules, ::Type{<:Number}) = rand() ≤ interactionrules.plysogeny
_lysogenic(phage, bact, interactionrules, ::Type{<:AbstractVector{Bool}}) = interactionrules.plysogeny[species(phage)]
_lysogenic(phage, bact, interactionrules, ::Type{<:AbstractVector}) = rand() ≤ interactionrules.plysogeny[species(phage)]

"""
lysogenic(phage::AbstractPhage, bact::AbstractBacterium,
                    interactionrules::AbstractInteractionRules)

Determines whether a phage will enter a lysogentic phase with the host as a
prophage or whether it will enter a lytic phase and kill its host immediately.
"""
function lysogenic(phage::AbstractPhage, bact::AbstractBacterium,
                    interactionrules::AbstractInteractionRules)
    return _lysogenic(phage, bact, interactionrules, plysogenytype(interactionrules))
end

"""
    lyses(bact::AbstractBacterium, grid::BactGrid, I::CartesianIndex,
                                interactionrules::AbstractInteractionRules)

Determines whether a bacterium with a latent phage will lyse. Is dependent of
the local density of the bacterial species and the `interactionrules`.
"""
function lyses(bact::AbstractBacterium, interactionrules::AbstractInteractionRules,
                                                                model)
    # compute the local density
    Rqs = interactionrules.Rqs
    # if radius is less than or equal to 1, lysis probability does not depend
    # on density
    return rand() < interactionrules.plysis
    Rqs ≤ 1 && return rand() < interactionrules.plysis
    sp = species(bact)
    ρ = density(bact, model, Rqs)
    # so the probability of lysis is proportional to the density
    # THE BACTERIUM ITSELF INCLUDED!
    return rand() < (ρ * interactionrules.plysis)
end

# checking whether a phage can infect a bacterium
infects(phage, bact, ir::InteractionRules, ::Type{Bool}) = ir.Pinf
infects(phage, bact, ir::InteractionRules, ::Type{<:Number}) = rand() ≤ ir.Pinf
infects(phage, bact, ir::InteractionRules, ::Type{<:AbstractMatrix{Bool}}) = ir.Pinf[species(bact), species(phage)]
infects(phage, bact, ir::InteractionRules, ::Type{<:AbstractMatrix{T}}) where {T<:Number} = rand() ≤ ir.Pinf[species(bact), species(phage)]

"""
    infects(phage, bact, ir::InteractionRules)

Test whether `phage` can infect `bact` according to the interaction rules `ir`.
Can be deterministic or stochastic, depending on the rules. Bacteria with a prophage
cannot be infected.
"""
function infects(phage::AbstractPhage, bact::AbstractBacterium,
                    ir::InteractionRules)
    haslatent(bact) && return false
    return infects(phage, bact, ir::InteractionRules, infectiontype(ir))
end

# default burstsize is Poisson distributed
# more general behaviour can also be implemented by extending these cases
burstsize(phage, bact, ir::InteractionRules, bstype) = throw(MethodError(burstsize, bact, ir, bstype))
burstsize(phage, bact, ir::InteractionRules, ::Type{<:Number}) = rand(Poisson(energy(bact) * ir.burstsize))

function burstsize(phage::AbstractPhage, bact::AbstractBacterium, ir::InteractionRules)
    return burstsize(phage, bact, ir, burssizetype(ir))
end

"""
    infect!(phage::AbstractPhage, bact::AbstractBacterium,
            interactionrules::AbstractInteractionRules, model)

Infect with a given probability. This function both determines whether
the infection takes place and, if it does, it updates the agents by removing
the host and adding new phage particles.
"""
function infect!(phage::AbstractPhage, bact::AbstractBacterium,
            interactionrules::AbstractInteractionRules, model)
    pos = bact.pos
    # decide to go lytic or not
    if lysogenic(phage, bact, interactionrules)
        prophage!(bact, species(phage))
        kill_agent!(phage, model)
    else
        sp = species(phage)
        kill_agent!(bact, model)
        kill_agent!(phage, model)
        for i in 1:burstsize(interactionrules)
            id = nextid(model)
            add_agent!(Phage(id, pos, sp), model)
        end
    end
end

function lyse!(bact::AbstractBacterium, model)
    sp = prophage(bact)
    pos = bact.pos
    kill_agent!(bact, model)
    for i in 1:burstsize(ir(model.properties))
        id = nextid(model)
        add_agent!(Phage(id, pos, sp), model)
    end
end

"""
    burstsize(interactionrules::AbstractInteractionRules)

Compute the burstsize of an infected bacterium, sampled from a Poisson
distribution.
"""
burstsize(interactionrules::AbstractInteractionRules) = rand(Poisson(interactionrules.burstsize))
