using Distributions
using FoodWebs
using OdinSon

"""
# TODO: shuffling really matters to the goodness of fit ... why?
I am not clear on why this matters to the statistics of the fit (the formula for ρ needs to
be over the entire distribution, as the current formula assumes equal marginal
distributions). This is because the formula's used assume that the marginal distributions
are equal (that is the distribution of the upper triangle and lower triangle are the same).
I wonder how far just making sure that ρ is calculated from the actual pair distribution
would work. (It works much better, though if the scaling is skewed, (lets say `f = 0.1`) the
distribution can become wonky again).

So at an extreme if the upper triangle was all the negative values and bottom all the
positive values, and `f = 0` then the eigenvalues would simply be the diagonal entries, but
from simulations smallish pertubations off of this lead to large changes in the eigenvalue
distribution, hence why the random matrix bound becomes more elongated, but the true
distribution becomes more spread out.

When one adds the parameters taken from the same off diagonal entries, but shuffled so that
the marginal distributions are made equal, the formula's become underestimates of the
stability of the system, but it does capture the shape better, and the two formula converge
as `f -> 1` as you would expect.
"""

nrow = 3
ncol = 3
fs = linspace(0, 1, nrow*ncol)
S = 500
cmats = zeros(S, S, length(fs))
for (i, f) in enumerate(fs)
    cmats[:, :, i] = random_predprey(may_network(S, 0.25), Uniform(0, 5), f, shuffle_prey=false)
end

figure()
for i in 1:size(cmats, 3)
    subplot(nrow, ncol, i)
    FoodWebs.plot_eigdist(cmats[:, :, i])
    FoodWebs.add_eig_ellipse(FoodWebs.shuffle_pairs(cmats[:, :, i]), stroke="#76a432")
end
