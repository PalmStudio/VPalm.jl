
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