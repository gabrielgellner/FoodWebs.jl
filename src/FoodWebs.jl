module FoodWebs
using Distributions

include("topology.jl")
export may_network
export cascade_network, generalized_cascade_network
export niche_network, generalized_niche_network

include("randcmat.jl")
export random_predprey

include("utils.jl")
export trophic_levels

include("circular_laws.jl")
export ellipsis_sample

end
