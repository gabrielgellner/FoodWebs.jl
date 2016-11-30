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
would work. (It works much better, though if the scaling is skewed, (lets say `f` = 0.1) the
distribution can become wonky again).

Even more interesting is that for low values of `f` we get a better fit with the naive
formula, whereas for higher values of `f` the empirical formula does better. This suggests
that something clear most be happening to the underlying parameters that I could try to
integrate.
"""

cmat1 = random_predprey(may_network(500, 0.25), Uniform(0, 5), 0.1, shuffle_prey=false)
cmat2 = random_predprey(may_network(500, 0.25), Uniform(0, 5), 0.4, shuffle_prey=false)
cmat3 = random_predprey(may_network(500, 0.25), Uniform(0, 5), 0.8, shuffle_prey=false)
cmat4 = random_predprey(may_network(500, 0.25), Uniform(0, 5), 1.0, shuffle_prey=false)

figure()
subplot(221)
FoodWebs.plot_eigdist(cmat1)
subplot(222)
FoodWebs.plot_eigdist(cmat2)
subplot(223)
FoodWebs.plot_eigdist(cmat3)
subplot(224)
FoodWebs.plot_eigdist(cmat4)
