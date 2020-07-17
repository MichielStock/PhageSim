module PhageSim

using Agents

abstract type AbstractRules end

export AbstractRules

include("bacteria.jl")
include("phages.jl")
include("interactions.jl")
include("simulate.jl")
include("utils.jl")


end # module
