module FoodWebs

using Distributions: Beta, Distribution, Uniform, logpdf, sample
using Graphs: DiGraph, Graph, add_edge!, gdistances, has_path, is_connected, strongly_connected_components,
    vertices, weakly_connected_components
using LinearAlgebra: Diagonal, I, LinearAlgebra, Symmetric, diag, diagind, diagm, eigmax, eigvals
using Random: AbstractRNG, Random
using SpecialFunctions: gamma
using Statistics: cov, mean, quantile, std
using TOML: TOML

include("topology.jl")
export PredationMatrix
export food_chain_network
export may_network
export cascade_network, generalized_cascade_network
export generalized_niche_network, minimum_potential_niche_network, nested_hierarchy_network,
    niche_network
export probabilistic_niche_network, probabilistic_niche_probabilities

include("randcmat.jl")
export expected_community_matrix, random_predprey

include("stability.jl")
export numerical_abscissa, spectral_abscissa
export λ_stability, ν_stability

include("utils.jl")
export connectance, trophic_incoherence, trophic_levels
export basal_species, disconnected_from_basal, top_species
export orient_two_way_links, remove_cannibalism

include("properties.jl")
export cannibals, intermediate_species, omnivores, species_in_loops
export generality, vulnerability, mean_max_similarity, trophic_similarity
export characteristic_path_length, structural_properties

include("feasible.jl")
export FeasibleCommunities, CommunityPrior, community_matrix, coordinates, diagnose, dimension,
    fluxes, growth_rates, halfspaces, interaction_matrix, isbounded, logdensity, structurally_singular
export draw_growth_rates, feasible_growth_rates, sample_communities

include("motifs.jl")
export motif_census, motif_ids, motif_name, motif_roles, motif_zscores, rewire_links

include("compartments.jl")
export compartment_homogeneity, compartments, modularity_zscore, species_roles, web_modularity

include("data.jl")
export FoodWeb, foodweb, foodweb_ids, read_foodweb
public foodweb_cache

include("biomass.jl")
export biomass, biomass_ratios, is_inverted_pyramid, is_pyramid, pyramid_slope, steady_state_biomass

include("layout.jl")
export foodwebplot, foodwebplot!, trophic_layout
public edge_crossings

include("elliptic_law.jl")
export eigenplot, eigenplot!, elliptic_law
export offdiag_pairs, offdiag_pairs_nz, shuffle_pairs
public EllipticLaw

end
