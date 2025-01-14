using FoodWebs
using Test

# Test set for trophic levels
linchain = [0 1 0; 0 0 1; 0 0 0]
@test trophic_levels(linchain) ≈ [1.0, 2.0, 3.0]

omnchain = [0 1 1; 0 0 1; 0 0 0]
trophic_levels(omnchain)
@test trophic_levels(omnchain) ≈ [1.0, 2.0, 2.5]

omn_relabel = [0 0 1; 1 0 1; 0 0 0]
@test trophic_levels(omn_relabel) ≈ [2.0, 1.0, 2.5]

cannibal = [1 1 0; 0 1 1; 0 0 1]
@test trophic_levels(cannibal) ≈ [1.0, 2.0, 3.0]

coupled_chain = [0 0 0 0 0; 0 0 1 0 0; 1 0 0 0 0; 1 0 0 0 0; 0 0 0 1 0]
@test trophic_levels(coupled_chain) ≈ [3, 1, 2, 2, 1]

coupled_chain_omn = [0 0 0 0 0; 1 0 1 0 0; 1 0 0 0 0; 1 0 0 0 0; 0 0 0 1 0]
@test trophic_levels(coupled_chain_omn) ≈ [8 / 3, 1, 2, 2, 1]

# Test set for connectance
@test isapprox(connectance(may_network(500, 0.25)), 0.25, rtol = 1e-2)
@test isapprox(connectance(may_network(500, 0.53)), 0.53, rtol = 1e-2)
@test isapprox(connectance(may_network(500, 0.72)), 0.72, rtol = 1e-2)

@test isapprox(connectance(cascade_network(500, 0.25)), 0.25, rtol = 1e-2)
@test isapprox(connectance(cascade_network(500, 0.57)), 0.57, rtol = 1e-2)
@test isapprox(connectance(cascade_network(500, 0.91)), 0.91, rtol = 1e-2)

@test isapprox(connectance(niche_network(500, 0.11)), 0.11, rtol = 1e-2)
@test isapprox(connectance(niche_network(501, 0.67)), 0.67, rtol = 1e-2)
@test isapprox(connectance(niche_network(502, 0.81)), 0.81, rtol = 1e-2)
connectance(niche_network(501, 0.67))

# Resilience Tests
A1 = [-1.0 0; 0 -10]
A2 = [-1.0 15; 0 -10]
@test λ_stability(A1) == λ_stability(A2)
@test ν_stability(A1) ≈ -1.0
@test isapprox(ν_stability(A2), 3.25, rtol = 1e-2)
