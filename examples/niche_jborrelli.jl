# # The Niche Model
# Example translated from [Jonathan J. Borrelli](https://jjborrelli.github.io/post/the-niche-model/)

# We need a few libraries for this one

using Distributions
using LightGraphs

# There are many models out there that have been proposed to generate food web network
# structures. One of the most popular is the niche model, first proposed by Williams and
# Martinez (2000) in [Simple Rules Yield Complex Food Webs](https://www.cs.dartmouth.edu/~rockmore/Williams2000Nature.pdf).
# This model has two parameters: the number of species (N), and the the fraction of realized
# links among species (connectance, C).
#
# Species are arranged along a single dimensional “niche” axis by assigning each of the `N`
# species a value drawn randomly from 0 to 1 (n_i). Each species also has a randomly determined
# feeding range (r_i), whose center (c_i) is a randomly drawn value between 0 and ni. The width
# of the feeding range is found by multiplying the niche value of species i by a value drawn
# from a beta distribution with an expected value of 2C. To get the appropriate beta
# distribution we use  β ( 1 , 1 2 C − 1 ) . For example, when  C = 0.2 the expected value of
# ```rand(Beta(1, ((1/(2*.2)) - 1)), 1)``` is 0.4.
#
# The n_vals() function below takes the number of species and connectance as arguments, and
# returns named tuple with N rows. Each row has the niche value, feeding center, and feedin
# range for one species.

function n_vals(S, C)
    niche = rand(S)
    r = rand(Beta(1, ((1 / (2 * C)) - 1)), S) .* niche
    ci = rand(Uniform(r / 2, niche), S)
    return (niche = niche, ci = ci, r = r)
end
