"""
    petiole_allometries(petiole_rachis_ratio_mean, petiole_rachis_ratio_sd, rachis_length, width_base, height_base, cpoint_width_intercept, cpoint_width_slope, cpoint_height_width_ratio, rng)

"""
function petiole_allometries(petiole_rachis_ratio_mean, petiole_rachis_ratio_sd, rachis_length, width_base, height_base, cpoint_width_intercept, cpoint_width_slope, cpoint_height_width_ratio, rng)
    petiole_rachis_length_ratio = mean_and_sd(petiole_rachis_ratio_mean, petiole_rachis_ratio_sd; rng=rng)
    petiole_length = petiole_rachis_length_ratio * rachis_length
    insertion_deviation_angle = normal_deviation_draw(5.0, rng)
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

- `rachis_length`: Length of rachis (cm)
- `cpoint_width_intercept`: Intercept of linear function
- `cpoint_width_slope`: Slope of linear function
"""
function width_at_cpoint(rachis_length, cpoint_width_intercept, cpoint_width_slope)
    return linear(rachis_length, cpoint_width_intercept, cpoint_width_slope)
end


function c_point_angle(leaf_rank, cpoint_decli_intercept, cpoint_decli_slope, cpoint_angle_SDP)
    angle = linear(leaf_rank, cpoint_decli_intercept, cpoint_decli_slope)
    angle += normal_deviation_draw(cpoint_angle_SDP)
    angle = abs(angle)
    return leaf_rank < 3 ? 0.5 * angle : angle
end

"""
    petiole_height(relative_position, cpoint_height, leaf_base_height)

Compute height profile along the petiole.

# Arguments

- `relative_position`: Position along the petiole (0-1)
- `cpoint_height`: Height of the leaf section at C point
- `leaf_base_height`: Height at the base of the leaf
"""
function petiole_height(relative_position, cpoint_height, leaf_base_height)
    return leaf_base_height - (leaf_base_height - cpoint_height) * sqrt(relative_position)
end

"""
    petiole_width(cpoint_width, relative_position, leaf_base_width)

Compute width profile along the petiole.

# Arguments

- `relative_position`: Position along the petiole (0-1) 
- `cpoint_width`: Width of the leaf at C point
- `leaf_base_width`: Width at base of leaf
"""
function petiole_width(relative_position, cpoint_width, leaf_base_width)
    return leaf_base_width - (leaf_base_width - cpoint_width) * relative_position^0.17
end