using FoodWebs
using Test
using Aqua
using ExplicitImports: test_explicit_imports

@testset "FoodWebs.jl" begin
    @testset "Aqua" begin
        Aqua.test_all(FoodWebs)
    end
    @testset "ExplicitImports" begin
        # `public` was added in Julia 1.11, so on older versions names like
        # `LinearAlgebra.checksquare` are not recognized as public
        test_explicit_imports(
            FoodWebs; all_qualified_accesses_are_public = VERSION >= v"1.11"
        )
    end
    include("topology.jl")
    include("utils.jl")
    include("randcmat.jl")
    include("resilience.jl")
    include("ellipsis_laws.jl")
end
