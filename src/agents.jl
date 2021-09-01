
export isbacterium, species, model_step!, agent_step!, Agent, infmodel, init_model, model_step_rand_phages

using Distributions: Poisson

using DrWatson: @dict

mutable struct Agent <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}

    isbact::Bool
    species::Int
    growthprog::Float64  # only for bact
end

# TODO: might work in 3D?


isbacterium(a::Agent) = a.isbact
species(a::Agent) = a.species

unitvector(ϕ) = reverse(sincos(ϕ))
rand_dir(Δ) = Δ .* (randn(), randn())

# does not do anything

model_step!(model) = nothing

function model_step_rand_phages(nphages=20)
    return function model_step!(model)
        for _ in 1:rand(Poisson(nphages))
            add_agent!(model, false, rand(1:model.nphagesp), 0.0)
        end
    end
end

function agent_step!(agent, model)
    # TODO: might add competition in reproduction?
    sp = agent.species
    pos = agent.pos
    l = model.l
    if isbacterium(agent)
        # BACTERIA
        # --------
        # die
        if rand(model.rng) < model.pdie
            model.nbacteria -= 1
            kill_agent!(agent, model)
            return 
        end

        # go over all neighbors, both to check if there is space to move

        dx, dy = rand_dir(model.Δbact)
        walk!(agent, (dx, dy), model)
        # now check if this place was free
        for neighbor in nearby_agents(agent, model, 2l, exact=true)
            if neighbor.id != agent.id && isbacterium(neighbor)
                walk!(agent, (-dx, -dy), model)  # move back...
                break
            end
        end

        # growth, based on a logistic model
        agent.growthprog += model.ΔE * rand(model.rng)

        # infection?
        for neighbor in nearby_agents(agent, model, l, exact=true)
            if !isbacterium(neighbor) && model.infection(agent.species, neighbor.species)
                # succesfull infection!
                for _ in 1:rand(Poisson(model.burstsize))
                    add_agent!(pos, model, false, neighbor.species, 0.0)
                end
                kill_agent!(agent, model)
                kill_agent!(neighbor, model)
                model.nbacteria -= 1
                return
            end
        end

        # decide to reproduce?
        if agent.growthprog ≥ 1.0
            orient = 2pi * rand(model.rng)
            daughter_pos = pos .+ 2.05l .* unitvector(orient)
            extent = model.space.extent
            daughter_pos = (daughter_pos .+ extent) .% extent
            for neighbor in nearby_agents(daughter_pos, model, 2l, exact=true)
                if isbacterium(neighbor)
                    return  # no place for reproduction!
                end
            end
            # add two daughter cells
            add_agent!(daughter_pos, model, true, sp, 0.0)
            agent.growthprog = 0.0
            model.nbacteria += 1
        end
    else
        # PHAGES
        # ------

        # decay
        if rand(model.rng) < model.pdecay
            kill_agent!(agent, model)
            return
        end

        # move
        dir = rand_dir(model.Δphage)
        walk!(agent, dir, model)

    end
    return
end

# infection models

"""
    infmodel(P[, p=1])

Generates a function `probinf` to simulate an infection,
i.e. `pinf(sp_bact, sp_phage) = rand() ≤ P[sp_bact, sp_phage] * p`
"""
function infmodel(P, p=1)
    @assert 0 < p ≤ 1
    return (sp_bact, sp_phage) -> P[sp_bact, sp_phage] > 0 ? rand() ≤ P[sp_bact, sp_phage] * p : false
end

function infmodel(P::Matrix{Bool}, p)
    @assert 0 < p ≤ 1
    return (sp_bact, sp_phage) -> P[sp_bact,sp_phage] ? rand() ≤ p : false
end

infmodel(P::Matrix{Bool}) = (sp_bact, sp_phage) -> P[sp_bact,sp_phage]

function infmodel(p::Real)
    @assert 0 ≤ p ≤ 1
    return (sp_bact, sp_phage) -> rand() ≤ p
end

function init_model(extent=(20, 20), spacing=min(extent...)/10; nbacteria=20, nbactsp=1, nphages=100,
    nphagesp=1,
    l=0.5,
    seed=24,
    ΔE = 0.2,
    Δbact = l,
    Δphage = Δbact,
    pdie = 0.01,
    pdecay=0.05,
    burstsize=10,
    infection=infmodel(0.5),
    periodic=true
    )

    # make model and space
    space = ContinuousSpace(extent, spacing; periodic)
    model = ABM(Agent, space,
                properties=@dict(
                    ΔE,
                    l,
                    Δbact,
                    Δphage,
                    nbactsp,
                    nphagesp,
                    pdie,
                    pdecay,
                    burstsize,
                    infection, nbacteria
                ),
                rng = MersenneTwister(seed))

    # populate
    maxiter = 100_000
    for _ in 1:maxiter
        pos = (rand() * extent[1], rand() * extent[2])
        !isempty(nearby_ids(pos, model, 2l, exact=true)) && continue
        sp = rand(1:nbactsp)
        add_agent!(pos, model, true, sp, rand())
        nagents(model) == nbacteria && break
    end
    if nagents(model) < nbacteria
        @warn "Could not place $nbacteria bacteria, model contains only $(nagents(model)) bacteria"
        model.nbacteria = nagents(model)
    end

    for _ in 1:nphages
        pos = (rand() * extent[1], rand() * extent[2])
        sp = rand(1:nphagesp)
        add_agent!(pos, model, false, sp, 0.0)
    end

    return model
end
