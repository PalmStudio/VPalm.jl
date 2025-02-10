


"""
    petiole(parent_node, index, scale, rachis_length, zenithal_insertion_angle, zenithal_cpoint_angle, parameters)

Make a leaf petiole.

# Arguments 

- `parent_node`: the parent node on which the petiole will be attached
- `index`: the MTG index of the petiole
- `scale`: the MTG scale of the petiole
- `rachis_length`: the rachis length, used to feed allometries to compute the petiole dimensions
- `zenithal_insertion_angle`: petiole insertion angle
- `zenithal_cpoint_angle`: angle at the C point (tip of the petiole, starting point of the rachis)
- `parameters`: a list of parameters as a `Dict{String}`:
    - "leaf_base_width": the base width of the petiole
    - "leaf_base_height": the base heigth of the petiole
    - "cpoint_width_intercept": petiole width at the c-point intercept for linear interpolation
    - "cpoint_width_slope": petiole width at the c-point slope for linear interpolation
    - "cpoint_height_width_ratio": height to width ratio at the C point
    - "petiole_rachis_ratio_mean": the average value of the ratio between rachis length and petiole length
    - "petiole_rachis_ratio_sd": its standard deviation
    - "rachis_nb_segments": the number of segments used to discretize the petiole
"""
function petiole(parent_node, index, scale, rachis_length, zenithal_insertion_angle, zenithal_cpoint_angle, parameters; rng=Random.MersenneTwister(1))
    petiole = Node(parent_node, NodeMTG("/", "Petiole", index, scale))
    compute_properties_petiole!(
        petiole,
        zenithal_insertion_angle, rachis_length, zenithal_cpoint_angle,
        parameters["leaf_base_width"], parameters["leaf_base_height"], parameters["cpoint_width_intercept"],
        parameters["cpoint_width_slope"], parameters["cpoint_height_width_ratio"],
        parameters["petiole_rachis_ratio_mean"],
        parameters["petiole_rachis_ratio_sd"], parameters["rachis_nb_segments"];
        rng=rng
    )

    for p in 1:parameters["rachis_nb_segments"]
        petiole_segment_node = Node(petiole, NodeMTG(p == 1 ? "/" : "<", "PetioleSegment", p, 6))
        compute_properties_petiole_section(petiole, petiole_segment_node, p, parameters["rachis_nb_segments"])
    end

    return petiole
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


function c_point_angle(leaf_rank, c_point_decli_intercept, c_point_decli_slope, c_point_angle_sdp)
    angle = linear(leaf_rank, c_point_decli_intercept, c_point_decli_slope)
    angle += normal_deviation_draw(c_point_angle_sdp)
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

"""
    generate_petiole_mesh(params)

Create full petiole mesh based on parameters.

# Arguments


"""
function generate_petiole_mesh(petiole_length, width_base, height_base, width_cpoint, height_cpoint,
    insertion_angle, zenithal_cpoint_angle, azimuthal_angle, nb_sections=15)

    petiole_sections = properties_petiole_section(
        petiole_length, width_base, height_base, width_cpoint, height_cpoint,
        insertion_angle, zenithal_cpoint_angle, azimuthal_angle, nb_sections
    )
    # Create segments
    segments = []
    #! Type this vector

    for p in 1:petiole_sections
        # Create elliptical cross section
        petiole_segment = elliptical_cylinder(p.width / 2.0, p.height / 2.0, p.length) |> Meshes.Rotate(RotXY(p.zenithal_angle_local, p.azimuthal_angle_local))

        push!(segments, petiole_segment)
    end

    # Merge segments into final mesh
    final_mesh = merge(segments...)

    return final_mesh
end