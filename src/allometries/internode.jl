"""
    Internode diameter model

Computes the diameter of an internode at a given rank.

# Arguments

- `internode_number`: The number of the internode.
- `nb_internodes`: The total number of internodes.
- `stem_diameter`: The diameter of the stem at the base.
- `stem_base_shrinkage`: The shrinkage coefficient at the stem base.
- `stem_top_shrinkage`: The shrinkage coefficient at the stem top.
- `leaves_in_sheath`: The number of leaves in the sheath.
"""
function internode_diameter(internode_number, nb_internodes, stem_diameter, stem_base_shrinkage, stem_top_shrinkage, leaves_in_sheath)
    # Shrink trunk base
    diameter = stem_diameter * (1 - exp(-stem_base_shrinkage * internode_number))
    # Shrink trunk top
    frond_rank = nb_internodes - internode_number - leaves_in_sheath
    reduction_factor = max(0, min(1, 1 - exp(-stem_top_shrinkage * frond_rank)))
    return diameter * reduction_factor
end


"""
    Internode length model

Computes the length of an internode at a given rank. The internode length is computed using a quadratic function.
The objective is to have a internodes that are short and growing for the first emitted leaves (before `nb_internodes_before_planting`), 
and then getting to a stable "constant" height, and at the end for the youngest leaves, having nodes currently growing (smaller).

# Arguments

- `internode_number`: The number of the internode.
- `stem_height`: The height of the stem.
- `nb_leaves`: The number of leaves.
- `nb_internodes`: The total number of internodes.
- `internode_rank_no_expansion`: The rank of the internode that will not expand.
- `nb_internodes_before_planting`: The number of internodes before planting.
- `internode_final_length`: The final length of the internode.
"""
function internode_length(internode_number, nb_internodes, stem_height, internode_rank_no_expansion, nb_internodes_before_planting, internode_min_height)

    # Computing the height of the internode considering all heights possible (between the first and the nb_internodes_before_planting internodes that are more and more tall, then constant, and then decreasing again from internode_rank_no_expansion until the youngest leaves)
    internode_heigth_final = (stem_height - internode_min_height * ((nb_internodes_before_planting * (nb_internodes_before_planting + 1)) / (2 * (1 - nb_internodes_before_planting)) + nb_internodes_before_planting - nb_internodes_before_planting / (1 - nb_internodes_before_planting) + (internode_rank_no_expansion + 1) / 2 + internode_rank_no_expansion - (nb_internodes + nb_internodes_before_planting) / internode_rank_no_expansion)) /
                             ((nb_internodes_before_planting * (nb_internodes_before_planting + 1)) / ((nb_internodes_before_planting - 1) * 2) + nb_internodes_before_planting / (1 - nb_internodes_before_planting) + nb_internodes - internode_rank_no_expansion - 2 - (internode_rank_no_expansion + 1) / 2 - (nb_internodes + nb_internodes_before_planting) / internode_rank_no_expansion)

    a = (internode_min_height - internode_heigth_final) / (1 - nb_internodes_before_planting)
    b = internode_min_height - a

    c = (internode_min_height - internode_heigth_final) / internode_rank_no_expansion
    # d = internode_min_height - c * (nb_internodes + nb_internodes_before_planting)
    # d = internode_heigth_final - (nb_internodes + nb_internodes_before_planting - internode_rank_no_expansion) * c
    d = internode_heigth_final - (nb_internodes - internode_rank_no_expansion) * c

    if internode_number < nb_internodes_before_planting
        length = a * internode_number + b
    elseif internode_number < (nb_internodes - internode_rank_no_expansion)
        length = internode_heigth_final
    else
        length = c * internode_number + d
    end

    return length
end

function internode_length_old(internode_number, nb_internodes, nb_leaves, stem_height, internode_rank_no_expansion, nb_internodes_before_planting, internode_final_length)
    rank = nb_leaves - internode_rank_no_expansion # Number of leaves not in expansion
    S = nb_internodes - nb_leaves # Number of internodes without leaves. Should be equal to 0 at adult stage
    nb_noexp = nb_internodes_before_planting + S + rank # Total number of internodes not in expansion
    nb_total = nb_internodes_before_planting + S + nb_leaves # Total number of internodes
    coeff = 2 * stem_height / (nb_internodes_before_planting * (2 * S + nb_internodes_before_planting + 1))
    l = coeff * nb_internodes_before_planting
    c = (l - internode_final_length) /
        (
        nb_noexp^2 -
        nb_total^2 -
        2 * nb_noexp * (rank - nb_leaves)
    )
    b = -2 * c * nb_noexp
    a = l - b * nb_noexp - c * nb_noexp^2
    if internode_number < nb_internodes_before_planting
        return coeff * internode_number
    elseif (internode_number >= nb_internodes_before_planting) & (internode_number < nb_noexp)
        return l
    else
        return a + b * internode_number + c * internode_number^2
    end
end


"""
    phyllotactic_angle(phyllotactic_angle_mean, phyllotactic_angle_sd; rng=Random.MersenneTwister(1234))

Computes the phyllotactic angle (°) using an average angle and a standard deviation (random draw from a normal distribution).

# Arguments

- `phyllotactic_angle_mean`: The average phyllotactic angle (°).
- `phyllotactic_angle_sd`: The standard deviation of the phyllotactic angle (°).

# Optional arguments

- `rng`: The random number generator.
"""
function phyllotactic_angle(phyllotactic_angle_mean, phyllotactic_angle_sd; rng=Random.MersenneTwister(1234))
    return mean_and_sd(phyllotactic_angle_mean, phyllotactic_angle_sd; rng=rng)
end
