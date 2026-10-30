"""
    spectral_abscissa(M)

The largest real part of the eigenvalues of `M`. The equilibrium of a system with Jacobian `M`
is locally stable if this is negative, and its magnitude is the asymptotic return rate after a
perturbation (resilience).
"""
spectral_abscissa(M::AbstractMatrix) = maximum(real, eigvals(M))

"""
    numerical_abscissa(M)

The largest eigenvalue of the symmetric part `(M + M') / 2` of `M`. This is the maximum
instantaneous growth rate of a perturbation, so the equilibrium is reactive (some
perturbations initially grow) if this is positive (Neubert & Caswell 1997). It bounds the
[`spectral_abscissa`](@ref) from above, with equality when `M` is normal.
"""
numerical_abscissa(M::AbstractMatrix) = eigmax(Symmetric((M + M') / 2))

"""
    λ_stability(M)

Alias for [`spectral_abscissa`](@ref).
"""
const λ_stability = spectral_abscissa

"""
    ν_stability(M)

Alias for [`numerical_abscissa`](@ref).
"""
const ν_stability = numerical_abscissa
