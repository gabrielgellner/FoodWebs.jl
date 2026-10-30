module FoodWebs

using Distributions: Beta, Distribution, Uniform, sample
using Graphs: DiGraph, is_connected
using LinearAlgebra: I, LinearAlgebra, diag, diagind, diagm, eigvals
using Random: AbstractRNG, Random
using Statistics: cov, mean

include("topology.jl")
export PredationMatrix
export food_chain_network
export may_network
export cascade_network, generalized_cascade_network
export generalized_niche_network, niche_network

include("randcmat.jl")
export random_predprey

include("resilience.jl")
export λ_stability, ν_stability

include("utils.jl")
export connectance, trophic_levels
export basal_species, disconnected_from_basal, top_species

include("ellipsis_laws.jl")
export ellipsis_sample
export offdiag_pairs, offdiag_pairs_nz, shuffle_pairs

end
