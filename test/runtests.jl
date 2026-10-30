using FoodWebs
using Test
using Aqua

@testset "FoodWebs.jl" begin
    @testset "Aqua" begin
        Aqua.test_all(FoodWebs)
    end
    include("topology.jl")
    include("utils.jl")
    include("randcmat.jl")
    include("resilience.jl")
    include("ellipsis_laws.jl")
end
