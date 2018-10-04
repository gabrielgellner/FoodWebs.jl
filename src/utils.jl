connectance(adj) = sum(adj) / size(adj, 2) ^ 2
connectance(adj::PredationMatrix) = 2 * sum(adj) / adj.S ^ 2

function basal_species(pmat::PredationMatrix)
    col_sums = sum(pmat, dims = 1)
    return [i for i in 1:length(col_sums) if col_sums[i] == 0]
end

function top_species(pmat::PredationMatrix)
    row_sums = sum(pmat, dims = 2)
    return [i for i in 1:length(row_sums) if row_sums[i] == 0]
end

function trophic_levels(adj)
    # the matrix comes in "predator" form but the matrix algebra requires it to be in
    # "prey" form, also needs to be float since we use transition probabilites
    pmat = collect(1.0 * adj')
    for i in 1:size(pmat, 1)
        # remove any cannibalism as it messes up the tp calc
        pmat[i, i] = 0
    end

    for i in 1:size(pmat, 1)
        total_prey = sum(pmat[i, :])
        if total_prey != 0
            pmat[i, :] ./= total_prey
        end
    end

    return pinv(I - pmat) * fill(1.0, size(pmat, 1))
end
