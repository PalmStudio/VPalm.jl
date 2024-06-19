
"""
    stem_bending(stem_bending_mean, stem_bending_sd; rng=Random.MersenneTwister(1234))

Computes the stem bending (°) using an average bending and a standard deviation (random draw from a normal distribution).

# Arguments

- `stem_bending_mean`: The average stem bending (°).
- `stem_bending_sd`: The standard deviation of the stem bending (°).

# Optional arguments

- `rng`: The random number generator.
"""
function stem_bending(stem_bending_mean, stem_bending_sd; rng=Random.MersenneTwister(1234))
    return stem_bending_mean + randn(rng) * stem_bending_sd
end

"""
    stem_height(nb_leaves_emitted, initial_stem_height, stem_height_coefficient, internode_length_at_maturity, stem_growth_start)

Computes the stem height (m) at a given number of leaves emitted.

# Arguments

- `nb_leaves_emitted`: The number of leaves emitted from planting.
- `initial_stem_height`: The initial stem height at planting (m).
- `stem_height_coefficient`: The coefficient of the exponential function.
- `internode_length_at_maturity`: The internode length when the plant is mature (m).
- `stem_growth_start`: The number of leaves emitted at which the stem starts to grow (m). This is because the stem does not grow at the same rate at the beginning of the plant's life,
because it first grows more in diameter than in height.
- `stem_height_variation`: The variation of the stem height (m) due to the random draw from a normal distribution.

# Optional arguments

- `rng`: The random number generator.

# Details

The stem height is computed using an exponential function for the first `stem_growth_start` leaves emitted, and then a linear function for the remaining leaves emitted.

Note that the stem height can also be subject to some variability using `stem_height_variation`, simulating natural variations that might occur in real-world scenarios, but
this variability will never make the stem height go below 30% of the intial computed height.
"""
function stem_height(nb_leaves_emitted, initial_stem_height, stem_height_coefficient, internode_length_at_maturity, stem_growth_start, stem_height_variation; rng=Random.MersenneTwister(1234))
    if nb_leaves_emitted <= stem_growth_start
        stem_height = initial_stem_height * exp(stem_height_coefficient * nb_leaves_emitted)
    else
        stem_height = internode_length_at_maturity * (nb_leaves_emitted - stem_growth_start) + initial_stem_height * exp(stem_height_coefficient * stem_growth_start)
    end

    # Add some variability to the stem_height, simulating natural variations that might occur in real-world scenarios:
    stem_height = max(0.3 * stem_height, stem_height + stem_height_variation * randn(rng))
    # Note that we use max(0.3 * stem_height,...) to ensure that the stem height is always at least 30% of the maximum height.

    return stem_height
end




"""
    stem_diameter(rachis_reference_length, stem_diameter_max, stem_diameter_slope, stem_diameter_inflection, stem_diameter_variation)

Computes the stem basis diameter (m)

# Arguments

- `rachis_reference_length`: The reference length (m) taken as the length of the rachis of the rank 1 leaf.
- `stem_diameter_max`: The maximum stem diameter (m).
- `stem_diameter_slope`: The slope used for the exponential function.
- `stem_diameter_inflection`: The inflection point used for the exponential function
- `stem_diameter_variation`: The variation of the stem diameter (m) due to the random draw from a normal distribution.

# Optional arguments

- `rng`: The random number generator.

# Details

The stem basis diameter is computed using an exponential function.

Note that the stem diameter can also be subject to some variability using `stem_diameter_variation`, simulating natural variations that might occur in real-world scenarios, but
this variability will never make the stem diameter go below 30% of the intial computed diameter.
Also, we remove 60% of the stem diameter (capped at 0.3m) to account for snags.
"""
function stem_diameter(rachis_reference_length, stem_diameter_max, stem_diameter_slope, stem_diameter_inflection, stem_diameter_variation;
    rng=Random.MersenneTwister(1234))
    stem_diameter = stem_diameter_max / (1.0 + exp( -4 * stem_diameter_slope * (rachis_reference_length - stem_diameter_inflection) ))

    # Add some variability to the stem_diameter, simulating natural variations that might occur in real-world scenarios:
    stem_diameter = max(0.3 * stem_diameter, stem_diameter + stem_diameter_variation * randn(rng))
    # Note that we use max(0.3 * stem_diameter,...) to ensure that the stem diameter is always at least 30% of the maximum diameter.
    stem_diameter -= min(0.3, 0.6 * stem_diameter)
    # Note that we remove extra diameter estimation due to snags

    return stem_diameter
end

"""
    internode_length(i, nb_internodes, stem_height)

Computes the internode length (m) for each internode of the stem.

# Arguments

- `i`: The index of the internode.
- `stem_height`: The stem height (m). Calculated using `stem_height()` function.
- `max_rank_expansion`: The rank from which internodes are no longer expanding.
- `apical_internode_length`: The length of the internode around the appical meristem (m).
- `nb_internodes_at_planting`: The number of internodes at planting (estimated).

# Details

"""
function internode_length(i, stem_height, nb_internodes, max_rank_expansion, apical_internode_length, nb_internodes_at_planting)
    rank = nb_internodes - max_rank_expansion
    coef = 2 * stem_height / (nb_internodes_at_planting * (nb_internodes_at_planting + 1))
    l =  nb_internodes_at_planting * coef
    c = (l - apical_internode_length) /
        ((nb_internodes_at_planting + rank)^2 - (nb_internodes_at_planting - nb_internodes)^2 - 2 * (nb_internodes_at_planting + rank) * (rank - nb_internodes))
    b = -2 * c * (nb_internodes_at_planting + rank)
    a = l - b * (nb_internodes_at_planting + rank) - c * (nb_internodes_at_planting + rank)^2
    if i < nb_internodes_at_planting
        return coef * i
    elseif i >= nb_internodes_at_planting && i < nb_internodes_at_planting + rank
        return l
    else
        return a + b * i + c * i^2
    end
end

"""
    internode_diameter(i, nb_internodes, stem_diameter, stem_shrinking_coefficient, nb_leaves_in_sheath)

Computes the internode diameter (m) for each internode of the stem.

# Arguments

- `i`: The index of the internode.
- `nb_internodes`: The number of internodes of the stem.
- `stem_diameter`: The stem diameter (m). Calculated using `stem_diameter()` function.
- `stem_base_shrinking_coefficient`: The coefficient of the exponential function used to shrink the base of the stem.
- `stem_top_shrinking_coefficient`: The coefficient of the exponential function used to shrink the top of the stem.
- `nb_leaves_in_sheath`: The number of leaves in the sheath.

# Details

The diameter of each internode of the stem is calculated from the base value of the stem diameter. We then apply two exponential functions to shrink the base and the top of the stem. The higher the values of the shrinking coefficients, the sharper the shrinkage.
"""

function internode_diameter(i, nb_internodes, stem_diameter, stem_base_shrinking_coefficient, stem_top_shrinking_coefficient, nb_leaves_in_sheath)
    # First shrink the base of the stem
    reduction_factor = 1 - exp(-i * stem_base_shrinking_coefficient)
    internode_diameter = stem_diameter * reduction_factor
    # Then shrink the top of the stem
    internode_rank = nb_internodes - i - nb_leaves_in_sheath
    reduction_factor = max(0,min(1, 1 - exp(-internode_rank * stem_top_shrinking_coefficient)))
    return internode_diameter * reduction_factor
end

"""
    phyllotactic_angle(frond_phylo_m, frond_phylo_sd; rng=Random.MersenneTwister(1234))

"""

function phyllotactic_angle(frond_phylo_m, frond_phylo_sd; rng=Random.MersenneTwister(1234))
    return frond_phylo_m + randn(rng) * frond_phylo_sd
end