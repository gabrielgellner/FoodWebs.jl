connectance(adj) = sum(adj) / size(adj, 2) ^ 2
connectance(adj::PredationMatrix) = 2 * sum(adj) / adj.S ^ 2

function trophic_levels(adj)
    # the matrix comes in "predator" form but the matrix algebra requires it to be in
    # "prey" form
    pmat = collect(1.0 * adj)
    for i in 1:size(pmat, 1)
        total_prey = sum(adj[i, :])
        if total_prey != 0
            pmat[i, :] ./= total_prey
        end
    end
    return pinv(I - pmat) * fill(1.0, size(pmat, 1))
end
