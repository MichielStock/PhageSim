#=
Created on 19/08/2021 09:45:55
Last update: 20/08/2021

@author: Michiel Stock
michielfmstock@gmail.com

A continuous bacteria-phage model
=#

using Agents, Random
using InteractiveDynamics, CairoMakie
using DrWatson

mutable struct Agent <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}

    isbact::Bool
    species::Int
    growthprog::Float64  # only for bact
end


isbacterium(a::Agent) = a.isbact

unitvector(ϕ) = reverse(sincos(ϕ))
rand_dir(Δ) = Δ .* (randn(), randn())



model_step!(model) = nothing

function agent_step!(agent, model)
    sp = agent.species
    pos = agent.pos
    if isbacterium(agent)
        # BACTERIA
        l = 1
        
        # die
        if rand(model.rng) < model.pdie
            model.nbacteria -= 1
            kill_agent!(agent, model)
            return 
        end

        # go over all neighbors, both to check if there is space to move

        dir = rand_dir(model.Δbact)
        walk!(agent, dir, model)
        # now check if this place was free
        for neighbor in nearby_agents(agent, model, 2l)
            if neighbor.id != agent.id && isbacterium(neighbor)
                walk!(agent, -1 .* dir, model)  # move back...
                break
            end
        end

        # growth, based on a logistic model
        agent.growthprog += model.ΔE * rand(model.rng)

        # infection?
        for neighbor in nearby_agents(agent, model, l)
            if !isbacterium(neighbor) && model.pinf(agent.species, neighbor.species) ≥ rand()
                # succesfull infection!
                for _ in 1:model.burstsize
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
            for neighbor in nearby_agents(daughter_pos, model, l)
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
        # decay
        if rand(model.rng) < model.pdecay
            kill_agent!(agent, model)
            return
        end

        # moving
        dir = rand_dir(model.Δphage)
        walk!(agent, dir, model)

    end
    return
end

pinf(sp_bact, sp_phage) = sp_bact == sp_phage ? 0.25 : 0.0


function init_model(spacesize=(20, 20);
    nbacteria=20, nbactsp=1, nphages=100, nphagesp=1,
    dt=0.25,
    l=0.5,
    seed=24,
    ΔE = 0.2,
    Δbact = l,
    Δphage = Δbact,
    pdie = 0.01,
    pdecay=0.05,
    burstsize=10,
    pinf=pinf
    )

    # make model and space
    space = ContinuousSpace(spacesize, 0.1, periodic=true,)
    model = ABM(Agent, space,
                properties=@dict(
                    ΔE,
                    l,
                    Δbact,
                    Δphage,
                    pdie,
                    pdecay,
                    burstsize,
                    pinf, nbacteria, K
                ),
                rng = MersenneTwister(seed))

    # populate
    while nagents(model) < nbacteria
        pos = (rand() * spacesize[1], rand() * spacesize[2])
        !isempty(nearby_ids(pos, model, 2l)) && continue
        sp = rand(1:nbactsp)
        add_agent!(pos, model, true, sp, rand())
    end

    for _ in 1:nphages
        pos = (rand() * spacesize[1], rand() * spacesize[2])
        sp = rand(1:nphagesp)
        add_agent!(pos, model, false, sp, 0.0)
    end

    return model
end

# should run until no more phages in the system

model = init_model((50, 50), nbactsp=3, nphagesp=3, nbacteria=500, nphages=100, burstsize=10, ΔE=0.2, pdie=0.01, l=0.1, K=1000.0, pdecay=0.1, Δphage=.5)

agentcolor(agent) = [:red, :blue, :green, :yellow, :black, :pink][agent.species]
agentsize(agent) = isbacterium(agent) ? 15.0 : 2.0
agentmarker(agent) = isbacterium(agent) ? :circle : :^

abm_video(
    "bact_test.mp4",
    model,
    agent_step!,
    model_step!;
    ac=agentcolor,
    as=agentsize,
    title = "Bacteria model",
    frames = 200,
    spf = 2,
    framerate = 5,
)


