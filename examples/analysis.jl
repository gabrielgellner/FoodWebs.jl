using FoodWebs
using Distributions
using PyPlot

function sample(nsample::Int, f::Float64)
    deigs = zeros(nsample)
    # interestingly |linspace| doesn't produce an array but instead an
    # iterator but I can use it like I would expect. I am not sure what the
    # performance trade off is
    rvals = linspace(0.0, 8.0, nsample)
    for i in 1:nsample
        web = may_network(100, 0.25)
        cmat = random_predprey(web, Uniform(rvals[i], rvals[i] + 1.0), f)
        deigs[i] = maximum(real(eigvals(cmat)))
    end
    hcat(rvals, deigs)
end

data = sample(500, 0.8)
fig1 = plot(data[:, 1], data[:, 2])
