"""
    compute_properties_petiole!(
        petiole_node,
        insertion_angle, rachis_length, zenithal_cpoint_angle,
        width_base, height_base, cpoint_width_intercept,
        cpoint_width_slope, cpoint_height_width_ratio,
        petiole_rachis_ratio_mean,
        petiole_rachis_ratio_sd, nb_sections;
        rng=Random.MersenneTwister(1)
    )

Compute the dimensional properties of a petiole.

# Arguments

- `petiole_node`: the MTG Node of the petiole
- `insertion_angle`: the angle of insertion of the petiole on the stem (°)
- `rachis_length`: the length of the rachis (m)
- `zenithal_cpoint_angle`: the zenithal angle of the C point of the petiole, *i.e.* the tip (°)
- `width_base`: the width of the petiole at its base (m)
- `height_base`: the height of the petiole at its base (m)
- `cpoint_width_intercept`: the intercept of the linear function for the width at the C point (m)
- `cpoint_width_slope`: the slope of the linear function for the width at the C point
- `cpoint_height_width_ratio`: the ratio of the height to width at the C point
- `petiole_rachis_ratio_mean`: the mean ratio of the petiole to rachis length
- `petiole_rachis_ratio_sd`: the standard deviation of the ratio of the petiole to rachis length
- `nb_sections`: the number of sections discretizing the petiole
- `rng=Random.MersenneTwister(1)`: the random number generator
"""
function compute_properties_petiole!(
    petiole_node,
    insertion_angle, rachis_length, zenithal_cpoint_angle,
    width_base, height_base, cpoint_width_intercept,
    cpoint_width_slope, cpoint_height_width_ratio,
    petiole_rachis_ratio_mean,
    petiole_rachis_ratio_sd, nb_sections;
    rng=Random.MersenneTwister(1)
)

    petiole = petiole_allometries(petiole_rachis_ratio_mean, petiole_rachis_ratio_sd, rachis_length, width_base, height_base, cpoint_width_intercept, cpoint_width_slope, cpoint_height_width_ratio, rng)

    petiole_node[:length] = petiole.length
    petiole_node[:azimuthal_angle] = petiole.azimuthal_angle
    petiole_node[:width_base] = petiole.width_base
    petiole_node[:height_base] = petiole.height_base
    petiole_node[:width_cpoint] = petiole.width_cpoint
    petiole_node[:height_cpoint] = petiole.height_cpoint
    petiole_node[:zenithal_insertion_angle] = insertion_angle
    petiole_node[:zenithal_cpoint_angle] = zenithal_cpoint_angle

    petiole_node[:section_length] = petiole.length / nb_sections
    petiole_node[:section_insertion_angle] = (zenithal_cpoint_angle - insertion_angle) / (nb_sections - 1)
    # Note: We use nb_sections - 1 because the first section already has an angle, the insertion_angle
    return nothing
end


"""
    compute_properties_petiole_section!(petiole_node, section_node, index, nb_sections)

Compute the dimension of a petiole section based on the dimensions of the petiole.

# Arguments

- `petiole_node`: the MTG Node of the petiole
- `section_node`: the MTG Node of the section to be computed
- `index`: the index of the section on the petiole, from 1 at the base to `nb_sections`.
- `nb_sections`: the number of sections discretizing the petiole

# Details

The `petiole_node` should have the following attributes:

- `width_base`: the width of the petiole at its base (m)
- `height_base`: the height of the petiole at its base (m)
- `width_cpoint`: the width of the petiole at the C point (m)
- `height_cpoint`: the height of the petiole at the C point (m)
- `section_length`: the length of the petiole sections (m)
- `insertion_angle`: the angle of insertion of the petiole on the stem (°)
- `section_insertion_angle`: the zenithal angle of insertion between the petioles sections (°)
- `azimuthal_angle`: the azimuthal angle at the insertion (°)
"""
function compute_properties_petiole_section!(petiole_node, section_node, index, nb_sections)
    petiole_section = properties_petiole_section(
        index, nb_sections, petiole_node.width_base, petiole_node.height_base,
        petiole_node.width_cpoint, petiole_node.height_cpoint, petiole_node.section_length,
        petiole_node.zenithal_insertion_angle, petiole_node.section_insertion_angle,
        petiole_node.azimuthal_angle
    )

    section_node.width = petiole_section.width
    section_node.height = petiole_section.height
    section_node.length = petiole_section.length
    section_node.zenithal_angle = petiole_section.zenithal_angle
    section_node.azimuthal_angle = petiole_section.azimuthal_angle
end

"""
    properties_petiole_section(
        index, nb_sections, width_base, height_base,
        width_cpoint, height_cpoint, petiole_section_length,
        petiole_insertion_angle, petiole_section_insertion_angle,
        azimuthal_angle
    )

Compute the properties of each section of the petiole.

# Arguments

- `index`: The index of the section within all sections (1-nb_sections)
- `nb_sections`: The number of sections discretizing the petiole
- `width_base`: Width of the petiole at its base (m)
- `heigth_base`: Height of the petiole at its base (m)
- `width_cpoint`: Width of the petiole at the C point (tip of the petiole, *i.e.* transition point to rachis, m)
- `height_cpoint`: Height at the C point (m)
- `petiole_section_length`: The length of the petiole sections (m)
- `petiole_insertion_angle`: Zenithal angle of insertion between the petiole and the stipe (local angle, relative to the stipe, °)
- `petiole_section_insertion_angle`: The zenithal angle of insertion between the petioles sections (°)
- `azimuthal_angle`: Azimuthal angle at the insertion (°)

# Returns 

A vector of dimensions for each section, given as a named tuple:

- width: width of the section (m)
- height: height of the section (m)
- length: length of the section (m)
- zenithal_angle_local: zenithal angle of the section (local angle, relative to the previous, °)
- azimuthal_angle_local: azimuthal angle of the section (local angle, relative to the previous, °)
"""
function properties_petiole_section(
    index, nb_sections, width_base, height_base,
    width_cpoint, height_cpoint, petiole_section_length,
    petiole_insertion_angle, petiole_section_insertion_angle,
    azimuthal_angle
)
    relative_position = index / nb_sections

    if index == 1
        zenithal_angle = petiole_insertion_angle
        section_width = width_base
        section_height = height_base
    else
        zenithal_angle = petiole_section_insertion_angle
        section_width = petiole_width(relative_position, width_base, width_cpoint)
        section_height = petiole_height(relative_position, height_base, height_cpoint)
    end

    if index == 2
        deviation_angle = azimuthal_angle
    else
        deviation_angle = 0.0
    end

    return (; width=section_width, height=section_height, length=petiole_section_length, zenithal_angle=zenithal_angle, azimuthal_angle=deviation_angle)
end