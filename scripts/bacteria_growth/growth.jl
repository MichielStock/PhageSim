#=
Created on 17/08/2021 14:49:58
Last update: Thursday 19/08/2021

@author: Michiel Stock
michielfmstock@gmail.com

More realistic growth model
=#

using Agents, LinearAlgebra
using InteractiveDynamics
using CairoMakie # choose plotting backend

mutable struct Bacterium <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}
    species::Int
    length::Float64
    orientation::Float64
    growthprog::Float64

    # node positions/forces
    p1::NTuple{2,Float64}
    p2::NTuple{2,Float64}
    f1::NTuple{2,Float64}
    f2::NTuple{2,Float64}
end

function Bacterium(id, pos, species, l, φ, g)
    a = Bacterium(id, pos, species, l, φ, g, (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0))
    update_nodes!(a)  # yields correct values for positions and forces
    return a
end

# this does not change anything, just updates the coordinates of the poles
function update_nodes!(a::Bacterium)
    offset = 0.5 * a.length .* unitvector(a.orientation)
    a.p1 = a.pos .+ offset
    a.p2 = a.pos .- offset
end

unitvector(φ) = reverse(sincos(φ))
cross2D(a, b) = a[1] * b[2] - a[2] * b[1]

project_in_space(pos, extent) = @. (pos + extent) % extent


function model_step!(model)
    for a in allagents(model)
        # TODO: implement a decay change
        if a.growthprog ≥ 1
            # When a cell has matured, it divides into two daughter cells on the
            # positions of its nodes.
            p1 = project_in_space(a.p1, model.space.extent)
            p2 = project_in_space(a.p2, model.space.extent)
            add_agent!(p1, Bacterium, model, 0.0, a.species, a.orientation, 0.0)
            add_agent!(p2, Bacterium, model, 0.0, a.species, a.orientation, 0.0)
            kill_agent!(a, model)
        else
            # The rest lengh of the internal spring grows with time. This causes
            # the nodes to physically separate.
            uv = unitvector(a.orientation)
            internalforce = model.hardness * (a.length - a.growthprog) .* uv
            a.f1 = -1 .* internalforce
            a.f2 = internalforce
        end
    end
    # Bacteria can interact with more than one other cell at the same time, therefore,
    # we need to specify the option `:all` in `interacting_pairs`
    # TODO: here we need to distinguish between bacteria and phages
    for (a1, a2) in interacting_pairs(model, 2.0, :all)
        interact!(a1, a2, model)
    end
end

function agent_step!(agent::Bacterium, model::ABM)
    fsym, compression, torque = transform_forces(agent)
    new_pos = agent.pos .+ model.dt * model.mobility .* fsym
    move_agent!(agent, new_pos, model)
    agent.length += model.dt * model.mobility .* compression
    agent.orientation += model.dt * model.mobility .* torque
    agent.growthprog += model.dt * model.growthrate  # TODO: depends on energy of system, species etc
    update_nodes!(agent)
    return agent.pos
end

function interact!(a1::Bacterium, a2::Bacterium, model)
    n11 = noderepulsion(a1.p1, a2.p1, model)
    n12 = noderepulsion(a1.p1, a2.p2, model)
    n21 = noderepulsion(a1.p2, a2.p1, model)
    n22 = noderepulsion(a1.p2, a2.p2, model)
    a1.f1 = @. a1.f1 + (n11 + n12)
    a1.f2 = @. a1.f2 + (n21 + n22)
    a2.f1 = @. a2.f1 - (n11 + n21)
    a2.f2 = @. a2.f2 - (n12 + n22)
end

function noderepulsion(p1::NTuple{2,Float64}, p2::NTuple{2,Float64}, model::ABM)
    delta = p1 .- p2
    distance = norm(delta)
    if distance ≤ 1.0
        uv = delta ./ distance
        return (model.hardness * (1.0 - distance)) .* uv
    end
    return (0.0, 0.0)
end

function transform_forces(agent::Bacterium)
    # symmetric forces (CM movement)
    fsym = agent.f1 .+ agent.f2
    # antisymmetric forces (compression, torque)
    fasym = agent.f1 .- agent.f2
    uv = unitvector(agent.orientation)
    compression = dot(uv, fasym)
    torque = 0.5 * cross2D(uv, fasym)
    return fsym, compression, torque
end

isbacterium(agent) = agent isa Bacterium

space = ContinuousSpace((50, 20), 0.2; periodic = true)
model = ABM(
    Bacterium,
    space,
    properties = Dict(:dt => 0.005, :hardness => 1e2, :mobility => 1.0, :growthrate => 0.1)
)

add_agent!((16.5, 8.0), Bacterium, model, 1, 0.0, 0.3, 0.0)
add_agent!((28.5, 12.0), Bacterium, model, 2, 0.0, 0.0, 0.0)
add_agent!((14.5, 3.0), Bacterium, model, 3, 0.0, 0.0, 0.0)

adata = [(isbacterium, count)]
mdata = []

#results, _ = run!(model, agent_step!, model_step!, 5000, adata=adata)


#plot(results.count_isbacterium)


function cassini_oval(agent)
    t = LinRange(0, 2π, 50)
    a = agent.growthprog
    b = 1
    m = @. 2 * sqrt((b^4 - a^4) + a^4 * cos(2 * t)^2) + 2 * a^2 * cos(2 * t)
    C = sqrt.(m / 2)

    x = C .* cos.(t)
    y = C .* sin.(t)

    uv = reverse(sincos(agent.orientation))
    θ = atan(uv[2], uv[1])
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]

    bacteria = R * permutedims([x y])
    coords = [Point2f0(x, y) for (x, y) in zip(bacteria[1, :], bacteria[2, :])]
    scale(Polygon(coords), 0.5)
end

bacteria_color(b) = CairoMakie.RGBf0(b.species * 3.14 % 1, 0.2, 0.2)



abm_video(
    "bacteria.mp4", model, agent_step!, model_step!;
    am = cassini_oval, ac = bacteria_color,
    spf = 50, framerate = 30, frames = 500,
    title = "Growing bacteria"
)
