#=
Created on Wednesday 26 March 2020
Last update: Friday 17 July 2020

@author: Michiel Stock
michielfmstock@gmail.com

General utilies for modelling.
=#

export bacteria1, bacteria2, bacteria3, phages1, phages2, phages3
export plothp

# SHOWING THE RULES
# -----------------

rulestring(x, level=0) = string(x)

function rulestring(rules::AbstractRules, level=0)
    rs = level > 0 ? "\n" : ""
    rs *= "\t"^level * string(typeof(rules)) * ":"
    fields = fieldnames(typeof(rules))
    for fn in fieldnames(typeof(rules))
        rs *= "\n" * "\t"^level * "- $fn : " * rulestring(getfield(rules, fn), level+1)
    end
    return rs
end

import Base: show
"""
Just an overloading for printing the various rules used in PhageSim.
"""
show(io::IO, rules::AbstractRules) = rulestring(rules) |> print

# COLLECTING DATA
# ---------------

bacteria(a::AbstractAgent, sp::Int) = bacteria(a) && species(a) == sp
bacteria1(a) = bacteria(a, 1)
bacteria2(a) = bacteria(a, 2)
bacteria3(a) = bacteria(a, 3)

phages(a::AbstractAgent, sp::Int) = phages(a) && species(a) == sp
phages1(a) = phages(a, 1)
phages2(a) = phages(a, 2)
phages3(a) = phages(a, 3)


using AgentsPlots


mshape(a::AbstractBacterium) = :circle
mshape(a::AbstractPhage) = :hex
as(a::AbstractPhage) = 1
as(a::AbstractBacterium) = 3energy(a)
offset(a::AbstractBacterium) = (0.0, 0.0)
offset(a::AbstractPhage) = (0.3randn(), 0.3randn())

function mcolor(a::Union{Bacterium,Phage})
    species(a) == 1 && return :red
    species(a) == 2 && return :blue
    species(a) == 3 && return :green
end

plothp(model) = plotabm(
    model;
    offset = offset,
    ac = mcolor,
    am = mshape,
    as = as,
    scheduler = by_type((Bacterium, Phage), false),
    grid = false,
    size = (800, 600),
    showaxis = false,
    aspect_ratio = :equal,
)
