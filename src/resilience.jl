λ_stability(M) = maximum(real.(eigvals(M)))
#NOTE: this is \nu_stability, which is pretty hard to read even with unicode
ν_stability(M) = λ_stability((M + M') / 2)
