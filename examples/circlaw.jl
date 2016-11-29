using Distributions
using FoodWebs
using PyPlot

"""
# TODO: shuffling really matters to the goodness of fit ... why?
I am not clear on why this matters to the statistics of the fit (the formula for ρ needs to
be over the entire distribution, as the current formula assumes equal marginal
distributions). This is because the formula's used assume that the marginal distributions
are equal (that is the distribution of the upper triangle and lower triangle are the same).
I wonder how far just making sure that ρ is calculated from the actual pair distribution
would work. (It works much better, though if the scaling is skewed, (lets say f = 0.1) the
distribution can become wonky again)
"""

cmat = random_predprey(may_network(500, 0.25), Uniform(0, 5), 0.5, shuffle_prey=false)
FoodWebs.plot_eigdist(cmat)
