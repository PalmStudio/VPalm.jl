
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

"""

    normal_random_draw(rng=Random.MersenneTwister(1234))

Draw a random value from a normal distribution with a mean of 0 and a standard deviation of 1.

# Optional arguments

- `rng`: The random number generator.
"""

function normal_random_draw(rng=Random.MersenneTwister(1234))
    x = 0.
    y = 0.
    z = 0.
    while true
        x = (randn(rng) * 2) - 1
        y = (randn(rng) * 2) - 1
        z = x^2 + y^2
        if z < 1
            break
        end
    end

    return x * sqrt(-2 * log(z) / z)

end


"""

    normal_deviation_draw(sd, rng=Random.MersenneTwister(1234))

Draw a random value from a normal distribution with a given standard deviation.

# Arguments

- `sd`: The standard deviation of the normal distribution.

# Optional arguments

- `rng`: The random number generator.
"""

function normal_deviation_draw(sd, rng=Random.MersenneTwister(1234))
    return sd * normal_random_draw(rng)
end


"""

    linear(x, intercept, slope)

Compute a linear function at given `x` value.

# Arguments

- `x`: The input value.
- `intercept`: The intercept of the linear function.
- `slope`: The slope of the linear function.
"""

function linear(x, intercept, slope)
    return intercept + slope * x
end