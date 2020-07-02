#=
Created on Wednesday 26 March 2020
Last update: Wednesday 01 July 2020

@author: Michiel Stock
michielfmstock@gmail.com

Functionality to simulate bacterium-phage interactions.
=#

export initialize_model, PhageSimRules
export br, pr, ir
using Agents

struct PhageSimRules{BR<:AbstractBacteriaRules,PR<:AbstractPhageRules,
                            IR<:AbstractInteractionRules}
    bacteriarules::BR
    phagerules::PR
    interactionrules::IR
end

function PhageSimRules(;
    prepr=0.3,
    pmove=0.3,
    pdie=0.1,
    pdecay=0.05,
    R=1,
    Pinf=0.1,
    burstsize=5,
    plysogeny=false,
    plysis=.1,
    Rqs=0
    )
    br = BacteriaRules(prepr, pmove, pdie)
    pr = PhageRules(pdecay, R)
    ir = InteractionRules(Pinf, burstsize, plysogeny, plysis, Rqs)
    return PhageSimRules(br, pr, ir)
end

br(rules::PhageSimRules) = rules.bacteriarules
pr(rules::PhageSimRules) = rules.phagerules
ir(rules::PhageSimRules) = rules.interactionrules

"""
    initialize_phagesim_model(;
        n_bacteria = 100,
        n_phages = 100,
        dims = (50, 50),
        n_bact_sp = 1,
        n_phage_sp = 1)

Initializes a model on a grid.
"""
function initialize_model(rules::PhageSimRules;
    n_bacteria = 100,
    n_phages = 100,
    dims = (50, 50),
    n_bact_sp = 1,
    n_phage_sp = 1,
)
    @assert n_bacteria ≤ prod(dims) "Too many bacteria for the grid"
    space = GridSpace(dims, moore=true)
    model = AgentBasedModel(Union{Bacterium,Phage}, space,
                            scheduler=by_type(true, true), warn=false,
                            properties=rules)
    id = 0
    for _ in 1:n_bacteria
        id += 1
        species = rand(1:n_bact_sp)
        pos = pos = (rand(1:dims[1]), rand(1:dims[2]))
        bact = Bacterium(id, pos, species)
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



#TODO: make pretty print for this function
