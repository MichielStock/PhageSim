#=
Created on Wednesday 26 March 2020
Last update: Friday 17 July 2020

@author: Michiel Stock
michielfmstock@gmail.com

Functionality to simulate bacterium-phage interactions.
=#

export initialize_model, PhageSimRules
export br, pr, ir
using Agents

struct PhageSimRules{BR<:AbstractBacteriaRules,PR<:AbstractPhageRules,
                            IR<:AbstractInteractionRules} <: AbstractRules
    bacteriarules::BR
    phagerules::PR
    interactionrules::IR
end


"""
PhageSimRules(;
    pmovebact=0.3,
    Ediv=4.0,
    Eupdate=ConstantEnergyUpdate(0.2, 8.0),
    pdecay=0.1,
    pmovephage=1.0,
    Pinf=0.1,
    burstsize=2,
    plysogeny=false,
    plysis=0.5,
    Rqs=0
    )

Construct an instance of `PhageSimRules`, itself composed of an instance of
`AbstractBacteriaRules`, `AbstractPhageRules`, and `AbstractInteractionRules`
(can be accesses using `br`, `pr` and `ir`, respectively).

The following parameters can be set (look in header for default values):
- `pmovebact`: probability that a bacteria will move, either a global value or a
                    vector with probabilites per species;
- `Ediv`: the energy level of the bacteria when they will attempt to divide,
                    value or vector;
- `Eupdate`: rules of type `AbstractEnergyupdate` that determine how the
                    gain or loose energy;
- `pdecay`: decay probability for the phages (same for all);
- `pmovephage`: probability that a phage will move in its turn;
- `Pinf`: the probability that a bacterium of a certain species can be infected
                    by a phage of a certain species, either a value (Float for
                    probability or a Boolean) or a matrix with the rows
                    representing bacteria species and the columns the phages;
- `burstsize`: the average number of phages that are generated **per unit of
                    energy in its host bacterium**;
- `plysogeny`: the probability that a phage species will enter a lysogenic
                    cycle, global value or vector of Boolean values or
                    probabilities;
- `plysis`: the probability that a bacterium with a prophage will lyse;
- `Rsq`: the radius for quorum sensing that is used to determine lysing probabity,
                    when set to greater than 1 the local bacterial density will
                    be multiplied with `plysis` to determine lysis.
"""
function PhageSimRules(;
    pmovebact=0.3,
    Ediv=4.0,
    Eupdate=ConstantEnergyUpdate(0.2, 8.0),
    pdecay=0.1,
    pmovephage=1.0,
    Pinf=0.1,
    burstsize=2,
    plysogeny=false,
    plysis=0.5,
    Rqs=0
    )
    br = BacteriaRules(pmovebact, Ediv, Eupdate)
    pr = PhageRules(pdecay, pmovephage)
    ir = InteractionRules(Pinf, burstsize, plysogeny, plysis, Rqs)
    return PhageSimRules(br, pr, ir)
end

br(rules::PhageSimRules) = rules.bacteriarules
br(bactrules::BacteriaRules) = bactrules
pr(rules::PhageSimRules) = rules.phagerules
ir(rules::PhageSimRules) = rules.interactionrules

"""
    initialize_phagesim_model(;
        n_bacteria = 100,
        n_phages = 100,
        dims = (50, 50),
        n_bact_sp = 1,
        n_phage_sp = 1,
        Einit=2.0)

Initializes a model on a grid.
"""
function initialize_model(rules::PhageSimRules;
    n_bacteria = 100,
    n_phages = 100,
    dims = (50, 50),
    n_bact_sp = 1,
    n_phage_sp = 1,
    Einit=1.0,
)
    @assert n_bacteria ≤ prod(dims) "Too many bacteria for the grid"
    space = GridSpace(dims, metric=:chebyshev)
    model = AgentBasedModel(Union{Bacterium,Phage}, space,
                            scheduler=by_type(true, true), warn=false,
                            properties=rules)
    id = 0
    for _ in 1:n_bacteria
        id += 1
        species = rand(1:n_bact_sp)
        pos = pos = (rand(1:dims[1]), rand(1:dims[2]))
        bact = Bacterium(id, pos, species, Einit)
        add_agent!(bact, model)
        move_agent_single!(bact, model)
    end
    for _ in 1:n_phages
        id += 1
        species = rand(1:n_phage_sp)
        pos = (rand(1:dims[1]), rand(1:dims[2]))
        add_agent!(Phage(id, pos, species), model)
    end
    return model
end
