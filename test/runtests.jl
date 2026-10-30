using FoodWebs
using Test
using Aqua
using ExplicitImports: test_explicit_imports

@testset "FoodWebs.jl" begin
    @testset "Aqua" begin
        Aqua.test_all(FoodWebs)
    end
    @testset "ExplicitImports" begin
        test_explicit_imports(FoodWebs)
    end
    include("topology.jl")
    include("utils.jl")
    include("randcmat.jl")
    include("resilience.jl")
    include("ellipsis_laws.jl")
end
