#=
Created on Friday 21 Feb 2020
Last update: Monday 29 June 2020

@author: Michiel Stock
michielfmstock@gmail.com

Model the interactions between bacteria and their phages
=#


export AbstractInteractionRules, InteractionRules

abstract type AbstractInteractionRules end

"""
Rules for the interactions between bacteria and phages.
"""
struct InteractionRules{TI,TB,TLG} <: AbstractInteractionRules
    Pinf::TI  # matrix to determine the chance of infection host-phage
    burstsize::TB  # average number of virons
    plysogeny::TLG  # probability of entering the lysogentic cycle
    plysis::Float64  # probability that an infected bacterium lyses
    Rqs::Int  # radius for quorum sensing
    function InteractionRules(Pinf, burstsize::Number,
            plysogeny::Union{BP,AbstractVector{BP}}=false,
            plysis=1.0, Rqs=0)  where {BP <: Union{Bool, AbstractFloat}}
        @assert all(0 .≤ plysogeny .≤ 1) "All values for `plysogeny` should be Booleans or probabilities"
        @assert !(plysogeny isa AbstractVector) || size(Pinf, 2) == length(plysogeny) "`plysogeny` has incorrect size"
        @assert 0 ≤ plysis ≤ 1 "`plysis` should be a valid probability"
        @assert Rqs ≥ 0 "radius for quorum sensing `R` should be nonzero"
        new{typeof(Pinf),typeof(burstsize),typeof(plysogeny)}(Pinf, burstsize, plysogeny,plysis)
    end
end

"""Get the type of `Pinf`, used for dispatch."""
infectiontype(ir::InteractionRules{TI,TB,TLG}) where {TI<:Any,TB<:Any,TLG<:Any} = TI

"""Get the type of `burstsize`, used for dispatch."""
burssizetype(ir::InteractionRules{TI,TB,TLG}) where {TI<:Any,TB<:Any,TLG<:Any} = TB

"""
    lysogenic(bact::AbstractBacterium, phagetype::Int,
                    interactionrules::AbstractInteractionRules)

Determines whether a phage will enter a lysogentic phase with the host as a
prophage or whether it will enter a lytic phase and kill its host immediately.
"""
function lysogenic(bact::AbstractBacterium, phage::AbstractPhage,
                    interactionrules::AbstractInteractionRules)
    plysogeny = interactionrules.plysogeny
    plysogeny isa Bool && return plysogeny
    plysogeny isa Number && return rand() ≤ plysogeny
    plysogeny isa AbstractVector{Bool} && return plysogeny[phagetype]
    plysogeny isa AbstractVector && return rand() ≤ plysogeny[phagetype]
    error("behaviour of `plysogeny` not defined (type is $(typeof(plysogeny)))")
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
    Rqs ≤ 1 && return return rand() < interactionrules.plysis
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
infects(phage, bact, ir::InteractionRules, ::Type{<:AbstractMatrix{T}}) where {T<:Number} = ir.Pinf[species(bact), species(phage)] ≤ rand()

"""
    infects(phage, bact, ir::InteractionRules)

Test whether `phage` can infect `bact` according to the interaction rules `ir`.
Can be deterministic or stochastic, depending on the rules.
"""
function infects(phage::AbstractPhage, bact::AbstractBacterium,
                    ir::InteractionRules)
    return infects(phage, bact, ir::InteractionRules, infectiontype(ir))
end

# default burstsize is Poisson distributed
# more general behaviour can also be implemented by extending these cases
# TODO: make that this is generated directly with Poisson type, likely faster
burstsize(phage, bact, ir::InteractionRules, bstype) = throw(MethodError(burstsize, bact, ir, bstype))
burstsize(phage, bact, ir::InteractionRules, ::Type{Number}) = rand(Poisson(ir.burstsize))

burstsize(phage::AbstractPhage, bact::AbstractBacterium, ir::InteractionRules) = burstsize(phage, bact, ir, burssizetype(ir))


"""
    infect!(phage::AbstractPhage, bact::AbstractBacterium,
            interactionrules::AbstractInteractionRules, model)

Infect with a given probability. This function both determine
"""
function infect!(phage::AbstractPhage, bact::AbstractBacterium,
            interactionrules::AbstractInteractionRules, model)
    pos = bact.pos
    kill_agent!(bact, model)
    kill_agent!(phage, model)
    sp = species(phage)
    for i in 1:burstsize(phage, bact, interactionrules)
        id = nextid(model)
        add_agent!(Phage(id, pos, species), model)
    end
end


end

"""
Compute the burstsize of an infected bacterium, sampled from a Poisson
distribution.
"""
burstsize(bact::AbstractBacterium, phagetype,
            interactionrules::AbstractInteractionRules) = rand(Poisson(interactionrules.burstsize))
