# Cross-validation of the niche and cascade models against EcologicalNetworksDynamics.jl (#19).
#
# EcologicalNetworksDynamics.jl (END) implements the niche model of Williams & Martinez (2000)
# and the cascade model of Cohen & Newman (1985) independently. Its generators (`niche_model`,
# `cascade_model`) return consumer × resource matrices without acceptance and draw from the
# global RNG. This script transposes them to prey × predator and applies the acceptance of
# FoodWebs.jl: a connected web with `abs(L - C S^2) <= max(tol C S^2, 1/2)`, `tol = 0.03`. For
# the species richness and connectance of each Williams & Martinez web it draws `N` accepted webs
# from each package and compares the distributions of the twelve `structural_properties` with a
# two-sample Kolmogorov–Smirnov test and the difference of means in standard errors.
#
# END is licensed AGPL-3.0, so it lives in this separate environment and the test suite does not
# depend on it. Run from the repository root (N defaults to 1000):
#
#     julia --project=validation validation/end_models.jl
#     N=200 julia --project=validation validation/end_models.jl

using FoodWebs: FoodWebs, cascade_network, connectance, foodweb, niche_network, structural_properties
using EcologicalNetworksDynamics: EcologicalNetworksDynamics
using Graphs: DiGraph, is_connected
using Printf: @sprintf
using Random: Random, Xoshiro
using Statistics: mean, std

const END = EcologicalNetworksDynamics.Internals
const N = parse(Int, get(ENV, "N", "1000"))
const TOL = 0.03

const WEBS = (
    "skipwith_pond", "little_rock_lake", "bridge_brook_lake", "chesapeake_bay", "ythan_estuary",
    "coachella_valley", "st_martin_island",
)

# The acceptance of FoodWebs.jl (`accept_web` in src/topology.jl)
function accepted(A, S, C)
    target = C * S^2
    return is_connected(DiGraph(A)) && abs(count(A) - target) <= max(TOL * target, 0.5)
end

# `N` accepted END webs as prey × predator matrices
function end_webs(model, S, C)
    webs = Matrix{Bool}[]
    while length(webs) < N
        A = permutedims(Matrix{Bool}(model(S, C)))
        accepted(A, S, C) && push!(webs, A)
    end
    return webs
end

# The two-sample Kolmogorov–Smirnov statistic and its asymptotic p-value
function ks_test(x, y)
    n, m = length(x), length(y)
    points = sort!(vcat(x, y))
    D = maximum(abs(count(<=(t), x) / n - count(<=(t), y) / m) for t in points)
    D == 0 && return D, 1.0
    λ = (sqrt(n * m / (n + m)) + 0.12 + 0.11 / sqrt(n * m / (n + m))) * D
    p = 2 * sum((-1)^(k - 1) * exp(-2 * k^2 * λ^2) for k in 1:100)
    return D, clamp(p, 0.0, 1.0)
end

function compare(model_name, ours, theirs)
    P, Q = structural_properties.(ours), structural_properties.(theirs)
    rows = []
    for name in keys(first(P))
        x = filter(!isnan, [p[name] for p in P])
        y = filter(!isnan, [q[name] for q in Q])
        D, p = ks_test(x, y)
        se = sqrt(var_or_zero(x) / length(x) + var_or_zero(y) / length(y))
        z = se == 0 ? 0.0 : (mean(x) - mean(y)) / se
        push!(
            rows, (;
                model = model_name, property = name, ours = mean(x), theirs = mean(y), z, D, p,
                nan = (length(P) - length(x), length(Q) - length(y)),
            )
        )
    end
    return rows
end
var_or_zero(x) = length(x) > 1 ? std(x)^2 : 0.0

Random.seed!(2026)
results = []
for (k, web) in enumerate(WEBS)
    fw = foodweb("williams_martinez_2000/$web")
    S, C = size(fw, 1), connectance(fw)
    rng = Xoshiro(k)
    for (name, ours, theirs) in (
            ("niche", () -> niche_network(rng, S, C; tol = TOL), END.niche_model),
            ("cascade", () -> cascade_network(rng, S, C; tol = TOL), END.cascade_model),
        )
        name == "cascade" && C > (S - 1) / (2S) && continue
        t = @elapsed rows = compare(name, [ours() for _ in 1:N], end_webs(theirs, S, C))
        for r in rows
            push!(results, (; web, S, C, r...))
        end
        @info "$web ($S species, C = $(round(C, digits = 3))), $name model: $(round(t, digits = 1)) s"
    end
end

tests = length(results)
println("$N webs per package and setting; $tests comparisons (property × web × model).\n")
println("| web | S | C | model | property | FoodWebs mean | END mean | z | KS D | KS p |")
println("|---|---|---|---|---|---|---|---|---|---|")
for r in sort(results; by = r -> r.p)[1:10]
    println(
        "| $(r.web) | $(r.S) | $(@sprintf("%.3f", r.C)) | $(r.model) | $(r.property) | ",
        "$(@sprintf("%.4g", r.ours)) | $(@sprintf("%.4g", r.theirs)) | $(@sprintf("%.2f", r.z)) | ",
        "$(@sprintf("%.3f", r.D)) | $(@sprintf("%.3g", r.p)) |"
    )
end
println()
println("KS p below 0.05: $(count(r -> r.p < 0.05, results)) of $tests (about $(round(Int, 0.05tests)) expected by chance)")
println("KS p below 0.01: $(count(r -> r.p < 0.01, results)) of $tests (about $(round(0.01tests, digits = 1)) expected)")
println("Smallest KS p: $(@sprintf("%.3g", minimum(r -> r.p, results))); Bonferroni threshold at 0.05: $(@sprintf("%.3g", 0.05 / tests))")
println("|z| above 3: $(count(r -> abs(r.z) > 3, results)) of $tests")
