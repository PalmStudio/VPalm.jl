
"""
    logistic(x, max, slope, inflection)

Compute a logistic function.

# Arguments

- `x`: The input value.
- `max`: The maximum value of the logistic function.
- `slope`: The slope of the logistic function.
- `inflection`: The inflection point of the logistic function.
"""

function logistic(x, max, slope, inflection)
    return max / (1. + exp(-4*slope * (x - inflection)))
end


"""

    mean_and_sd(mean, sd; rng=Random.MersenneTwister(1234))

Compute a random value from a normal distribution with a given mean and standard deviation.

# Arguments

- `mean`: The mean of the normal distribution.
- `sd`: The standard deviation of the normal distribution.

# Optional arguments

- `rng`: The random number generator.
"""

function mean_and_sd(mean, sd; rng=Random.MersenneTwister(1234))
    return mean + randn(rng) * sd
end