"""
    petiole_allometries(petiole_rachis_ratio_mean, petiole_rachis_ratio_sd, rachis_length, width_base, height_base, cpoint_width_intercept, cpoint_width_slope, cpoint_height_width_ratio, rng)

Compute the petiole allometries.

# Arguments

- `petiole_rachis_ratio_mean`: Average value of the petiole/rachis length ratio
- `petiole_rachis_ratio_sd`: Standard deviation of the petiole/rachis length ratio
- `rachis_length`: Length of the rachis (m)
- `width_base`: Width at the base of the leaf (m)
- `height_base`: Height at the base of the leaf (m)
- `cpoint_width_intercept`: Intercept of the linear relationship between rachis width at C point and rachis length (m)
- `cpoint_width_slope`: Slope of the linear relationship
- `cpoint_height_width_ratio`: Ratio between the height and width of the leaf at C point
- `rng`: Random number generator

# Returns

A named tuple with the following keys:

- `length`: Length of the petiole (m)
- `azimuthal_angle`: Azimuthal angle of the petiole (°)
- `width_base`: Width at the base of the petiole (m)
- `height_base`: Height at the base of the petiole (m)
- `width_cpoint`: Width at the C point of the petiole (m)
- `height_cpoint`: Height at the C point of the petiole (m)

# Details

The petiole length is computed as the product of the rachis length and the petiole/rachis length ratio. 
The azimuthal angle is drawn from a normal distribution with a standard deviation of 5°. 
The width at the C point is computed using a linear function of the rachis length. 
The height at the C point is computed as the product of the width at the C point and the height/width ratio.
"""
function petiole_allometries(petiole_rachis_ratio_mean, petiole_rachis_ratio_sd, rachis_length, width_base, height_base, cpoint_width_intercept, cpoint_width_slope, cpoint_height_width_ratio, rng)
    petiole_rachis_length_ratio = mean_and_sd(petiole_rachis_ratio_mean, petiole_rachis_ratio_sd; rng=rng)
    petiole_length = petiole_rachis_length_ratio * rachis_length
    insertion_deviation_angle = normal_deviation_draw(5.0u"°", rng) #! this should be a parameter. And we should be able to remove the randomness with an option.
    width_cpoint = width_at_cpoint(rachis_length, cpoint_width_intercept, cpoint_width_slope)
    height_cpoint = cpoint_height_width_ratio * width_cpoint
    #! These should be allometries relative to leaf length, because tiny leaves don't have big bases:
    # width_base = width_base 
    # height_base = height_base
    return (length=petiole_length, azimuthal_angle=insertion_deviation_angle, width_base=width_base, height_base=height_base, width_cpoint=width_cpoint, height_cpoint=height_cpoint)
end

"""
    width_at_cpoint(rachis_length, cpoint_width_intercept, cpoint_width_slope)

Compute width at C point based on rachis length.

# Arguments

- `rachis_length`: Length of rachis (m)
- `cpoint_width_intercept`: Intercept of linear function (m)
- `cpoint_width_slope`: Slope of linear function
"""
function width_at_cpoint(rachis_length, cpoint_width_intercept, cpoint_width_slope)
    return linear(rachis_length, cpoint_width_intercept, cpoint_width_slope)
end


"""
    c_point_angle(leaf_rank, cpoint_decli_intercept, cpoint_decli_slope, cpoint_angle_SDP; rng)

Compute the angle at the C point of the leaf.

# Arguments

- `leaf_rank`: Rank of the leaf
- `cpoint_decli_intercept`: Intercept of the linear relationship between leaf rank and C point declination
- `cpoint_decli_slope`: Slope of the linear relationship
- `cpoint_angle_SDP`: Standard deviation of the C point angle
- `rng`: Random number generator

# Returns

The zenithal angle at the C point of the leaf (°)
"""
function c_point_angle(leaf_rank, cpoint_decli_intercept, cpoint_decli_slope, cpoint_angle_SDP; rng=Random.MersenneTwister(1))
    angle = linear(leaf_rank, cpoint_decli_intercept, cpoint_decli_slope)
    angle += normal_deviation_draw(cpoint_angle_SDP, rng) |> abs
    angle = leaf_rank < 3 ? 0.5 * angle : angle
    return angle * unit(cpoint_decli_intercept)
end

"""
    petiole_height(relative_position, height_cpoint, height_base)

Compute the height profile along the petiole (m).

# Arguments

- `relative_position`: Position along the petiole (0-1)
- `height_base`: Height at the base of the leaf
- `height_cpoint`: Height of the leaf section at C point
"""
function petiole_height(relative_position, height_base, height_cpoint)
    return height_base - (height_base - height_cpoint) * sqrt(relative_position)
end

"""
    petiole_width(relative_position, width_cpoint, width_base)

Compute the width profile along the petiole (m).

# Arguments

- `relative_position`: Position along the petiole (0-1) 
- `width_base`: Width at base of leaf
- `width_cpoint`: Width of the leaf at C point
"""
function petiole_width(relative_position, width_base, width_cpoint)
    return width_base - (width_base - width_cpoint) * relative_position^0.17
end