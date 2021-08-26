#=
Created on
Last update:

@author: Michiel Stock
michielfmstock@gmail.com

Simple test of Agents.jl for host-phage interaction.
=#

using Agents, AgentsPlots, StatsPlots

mutable struct Bacterium <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    species::Int
end

mutable struct Phage <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    species::Int
end

function initialize_model(;
    n_bacteria = 100,
    n_phages = 200,
    dims = (50, 50),
    n_bact_sp = 1,
    n_phage_sp = 1,
)
    space = GridSpace(dims, moore = true)
    model =
        AgentBasedModel(Union{Bacterium,Phage}, space, scheduler = by_type(true, true), warn = false)
    id = 0
    for _ in 1:n_bacteria
        id += 1
        species = rand(1:n_bact_sp)
        bact = Bacterium(id, (0, 0), species)
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

function agent_step!(bact::Bacterium, model)
    prepr, pmove, pdie = 0.05, 0.2, 0.01
    rand() < pmove && move!(bact, model)
    rand() < prepr && reproduce!(bact, model)
    rand() < pdie && kill_agent!(bact, model)
end

function agent_step!(phage::Phage, model)
    pmove, pdie, pinfect = 0.75, 0.05, 0.5
    burstsize = 5
    rand() < pmove && move!(phage, model)
    agents = get_node_agents(phage.pos, model)
    hosts = filter!(x -> isa(x, Bacterium), agents)
    for host in hosts
        if phage.species == host.species && rand() < pinfect
            for _ in 1:burstsize
                add_agent!(Phage(nextid(model), phage.pos, phage.species), model)
            end
            kill_agent!(host, model)
            kill_agent!(phage, model)
            return
        end
    end
    rand() < pdie && kill_agent!(phage, model)
    return
end

function move!(bact::Bacterium, model)
    neighbors = node_neighbors(bact, model)
    node = rand(neighbors)
    if isfree(node, model)
        move_agent!(bact, node, model)
    end
end

function move!(phage::Phage, model)
    r = 2  # phage can move two positions
    neighbors = node_neighbors(phage, model)
    node = rand(neighbors)
    move_agent!(phage, node, model)
end

isfree(node, model) = count(bacteria, get_node_agents(node, model)) == 0

function reproduce!(bact::Bacterium, model)
    neighbors = node_neighbors(bact, model)
    node = rand(neighbors)
    if isfree(node, model)
        id = nextid(model)
        add_agent!(Bacterium(id, node, bact.species), model)
    end
end



mshape(a::Bacterium) = :circle
mshape(a::Phage) = :hex
as(a::Phage) = 2
as(a::Bacterium) = 10
offset(a::Phage) = (0.1randn(), 0.1randn())
offset(a::Bacterium) = (0.0, 0.0)

function mcolor(a::Union{Bacterium,Phage})
    a.species == 1 && return :red
    a.species == 2 && return :blue
    a.species == 3 && return :green
end

plotabm(
    model;
    offset = offset,
    ac=mcolor,
    am = mshape,
    as = as,
    scheduler = by_type((Bacterium, Phage), false),
    grid = false,
    size = (800, 600),
    showaxis = false,
    aspect_ratio = :equal,
)

bacteria(a) = typeof(a) == Bacterium
phages(a) = typeof(a) == Phage
adata = [(bacteria, count), (phages, count)]

n_steps = 1000
model = initialize_model(n_bact_sp=3, n_phage_sp=3)

results, _ = run!(model, agent_step!, n_steps, adata = adata)

plot(results.step, results.count_bacteria, label=:bacteria)
plot!(results.step, results.count_phages, label=:phages)

model = initialize_model(n_bact_sp=3, n_phage_sp=3, dims=(20, 20), n_bacteria=30, n_phages=40)
anim = @animate for i in 0:250
    i > 0 && step!(model, agent_step!)
    p1 = plotabm(
        model;
        offset = offset,
        ac=mcolor,
        am = mshape,
        as = as,
        scheduler = by_type((Bacterium, Phage), false),
        grid = false,
        size = (800, 600),
        showaxis = false,
        aspect_ratio = :equal,
    )
    title!(p1, "step $(i)")
end

gif(anim, "bactphage.gif", fps = 10)
