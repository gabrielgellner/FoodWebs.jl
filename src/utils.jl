
function trophic_levels(adj::Matrix)
    # only look at the upper triangle to avoid double counting
    uadj = triu(adj)
    S = size(adj, 1)
    @assert S == size(uadj, 2) # only square matrices
    maxlevels = 8 # assume no trophic positions higher than 8 TODO: remove the need for this
    # basal species start with tl 0 all other levels are left undetermined
    NA = -100 # sentinal to mean not defined.
    tl = Vector{Int}(S)
    for (i, sp) = enumerate(1:S)
        if sum(uadj[:, sp]) == 0
            tl[i] = 0
        else
            tl[i] = NA
        end
    end

    # assign integer trophic levels
    for sp = 1:S, level = 1:maxlevels
        for resource in find(tl .== (level - 1))
            if uadj[resource, sp] == 1 && uadj[sp, resource] == 0
                tl[sp] = level
            end
        end
    end

    # adjust for omnivory
    bylevel = Matrix{Int}(S, S)
    for i = 1:S, j = 1:S
        bylevel[i, j] = uadj[i, j]*(tl[i] + 1) + 1
        if bylevel[i, j] == 1
            bylevel[i, j] = 0
        end
    end

    # calculate the actual trophic position
    tp = zeros(S)
    for sp = 1:S
        if sum(uadj[:, sp]) == 0 # basal species
            tp[sp] = 1
        else
            # assume that all feeding interations are equal so just sum the number of
            # interactions for each species and divide to get the binary trophic position
            tp[sp] = sum(bylevel[:, sp])/sum(uadj[:, sp])
        end
    end
    return tp
end
