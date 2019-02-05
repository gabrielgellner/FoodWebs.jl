module FoodWebs

using LinearAlgebra
using Distributions
using LightGraphs

include("topology.jl")
export food_chain_network
export may_network
export cascade_network, generalized_cascade_network
export niche_network, generalized_niche_network

include("randcmat.jl")
export random_predprey

include("resilience.jl")
export λ_stability, ν_stability

include("utils.jl")
export trophic_levels, connectance
export basal_species, top_species

include("ellipsis_laws.jl")
export ellipsis_sample

end
