using Distributions
using FoodWebs
using PyPlot
#TODO: shuffling really matters to the goodness of fit ... why? I am not clear on why this matters
# to the statistics of the fit
cmat = random_predprey(may_network(500, 0.25), Uniform(0, 5), 0.1, shuffle_prey=false)
elaw = ellipsis_sample(cmat)
mean(diag(cmat))
elaw.μ_diag
FoodWebs.plot_eigs(cmat)
FoodWebs.eig_ellipse(cmat)
