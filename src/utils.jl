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
    return max / (1. + exp(-4 * slope * (x - inflection)))
end

function logistic(x::T, max, slope, inflection::T) where T<:Quantity
    return max / (1. + exp(-4 * slope * (ustrip(x) - ustrip(inflection))))
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

    normal_deviation_draw(sd, rng=Random.MersenneTwister(1234))

Draw a random value from a normal distribution with a given standard deviation.

# Arguments

- `sd`: The standard deviation of the normal distribution.

# Optional arguments

- `rng`: The random number generator.
"""
function normal_deviation_draw(sd, rng=Random.MersenneTwister(1234))
    return sd * rand(rng)
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


"""
    exponetial(x, a, b)

Compute an exponential function at given `x` value.

# Arguments

- `x`: The input value.
- `a`: The coefficient `a` of the exponential function.
- `b`: The coefficient `b` of the exponential function.

# Note

The exponential function is defined as `a * exp(b * x)`.
"""
function exponetial(x, a, b)
    return a * exp(b * x)
end

"""
    normal_deviation_percent_draw(value, sdp, rng)

Calculate a normally distributed random deviation based on a percentage of the value.

# Arguments
- `value`: Base value.
- `sd`: Standard deviation in %.
- `rng`: Random number generator.

# Returns

- The random deviation.
"""
function normal_deviation_percent_draw(value, sd, rng)
    return normal_deviation_draw(sd, rng) * (value / 100.0)
end

"""
    beta_distribution_norm(x, xm, ym)

Calculate the normalized beta distribution value at point x.

# Arguments
- `x`: Position [0 to 1].
- `xm`: Mode of the beta distribution.
- `ym`: Maximum value of the beta distribution.

# Returns
- Normalized beta distribution value.
"""
function beta_distribution_norm(x, xm, ym)
    if x <= 0 || x >= 1
        return 0.0
    end

    # Bell curve approximation
    return ym * (x / xm) * ((1 - x) / (1 - xm)) * 4
end