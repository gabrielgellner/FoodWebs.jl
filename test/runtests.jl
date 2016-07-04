using FoodWebs
using Base.Test

#TODO: upgrade to testsets once 0.5 is out
chain = [0 1 0; 1 0 1; 0 1 0]
@test_approx_eq trophic_levels(chain) [1.0 2.0 3.0]
omnchain = [0 1 1; 1 0 1; 1 1 0]
@test_approx_eq trophic_levels(omnchain) [1.0, 2.0, 2.5]
