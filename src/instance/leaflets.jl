mutable struct Leaflet
    group::Int
    group_size::Int
    plane::Int
end

function leaflets(rachis_node, index, scale, leaf_rank, rachis_length, height_cpoint, width_cpoint, zenithal_cpoint_angle, parameters; rng=rng)
    nb_leaflets = compute_number_of_leaflets(rachis_length, parameters["leaflets_nb_max"], parameters["leaflets_nb_min"], parameters["leaflets_nb_slope"], parameters["leaflets_nb_inflexion"], parameters["nbLeaflets_SDP"]; rng=rng)

    leaflets_type_frequency = compute_leaflet_type_frequencies(parameters["leaflet_frequency_high"], parameters["leaflet_frequency_low"])

    # Structure of Arrays approach
    leaflets = (
        group=zeros(Int, nb_leaflets),
        group_size=zeros(Int, nb_leaflets),
        plane=zeros(Int, nb_leaflets)
    )

    group_leaflets!(leaflets, nb_leaflets, leaflets_type_frequency, parameters["leaflet_position_shape_coefficient"], rng)

    return leaflets
end

function group_leaflets!(leaflets, nb_leaflets, leaflets_type_frequency, shape_coefficient, rng)
    group = 1
    l = 1

    while l <= nb_leaflets
        relative_rank = l / nb_leaflets
        relative_position = relative_leaflet_position(relative_rank, shape_coefficient)

        group_size = draw_group_size(relative_position, leaflets_type_frequency, rng)
        group_size = min(group_size, nb_leaflets - l + 1)

        segment = clamp(floor(Int, relative_position * length(leaflets_type_frequency)), 1, length(leaflets_type_frequency))
        frequencies = leaflets_type_frequency[segment]

        # The first leaflet in the group is always plane=1
        leaflets.group[l] = group
        leaflets.group_size[l] = group_size
        leaflets.plane[l] = 1

        # Process other leaflets in the group if any
        for f in 1:(group_size-1)
            idx = l + f
            leaflets.group[idx] = group
            leaflets.group_size[idx] = group_size

            # Determine plane type based on frequencies
            if rand(rng) > (frequencies.medium / (frequencies.medium + frequencies.low))
                leaflets.plane[idx] = -1
            else
                leaflets.plane[idx] = 0
            end
        end

        l += group_size
        group += 1
    end

    return leaflets
end


"""
    compute_number_of_leaflets(rachis_final_length, nb_max, nb_slope, nb_infl, nbLeaflets_SDP; rng)

Compute the number of leaflets based on the logistic function, a standard deviation and a minimum value allowed.

# Arguments

- `rachis_final_length`: Final length of the rachis (m).
- `nb_max`: Maximum number of leaflets.
- `nb_min`: Minimum number of leaflets.
- `nb_slope`: Slope parameter for the logistic function (leaflet m⁻¹).
- `nb_infl`: Inflection point parameter for the logistic function (m).
- `nbLeaflets_SDP`: Standard deviation of the normal distribution for the number of leaflets.
- `rng`: Random number generator.

# Returns

The computed number of leaflets (integer).
"""
function compute_number_of_leaflets(rachis_final_length, nb_max, nb_min, nb_slope, nb_infl, nbLeaflets_SDP; rng)

    @assert rachis_final_length >= 0 "Rachis length must be non-negative"
    @assert nb_max > 0 "Maximum number of leaflets must be positive"
    @assert nb_min >= 0 "Minimum number of leaflets must be non-negative"
    @assert nb_slope > 0 "Slope parameter must be positive"
    @assert nb_infl > 0 "Inflection point must be positive"

    nb_leaflets = logistic(rachis_final_length, nb_max, nb_slope, nb_infl)

    deviation_factor = rachis_final_length < 1.0 ? 0.3 : 1.0
    nb_leaflets += deviation_factor * normal_deviation_draw(nbLeaflets_SDP, rng)

    nb_leaflets = round(Int, max(nb_min, nb_leaflets))

    return nb_leaflets
end

"""
    relative_leaflet_position(relative_rank, coef_dispo)

Compute the relative leaflet position on the rachis.

# Arguments

- `relative_rank`: Relative leaflet rank, usually in the form of (0 to 1].
- `shape_coefficient`: Shape coefficient (around 0).

# Returns

The relative leaflet position, in the same form as `relative_rank`, usually (0 to 1].
"""
function relative_leaflet_position(relative_rank, shape_coefficient)
    return ((1.0 + shape_coefficient) * (relative_rank^2)) / (1.0 + shape_coefficient * (relative_rank^2))
end


"""
    compute_leaflet_type_frequencies(leaflet_frequency_high, leaflet_frequency_low)

Compute the frequency of leaflet type within the sub-sections of a rachis.

# Arguments

- `leaflet_frequency_high`: Vector of frequency values for the +1 leaflet types (high) along the rachis sub-sections.
- `leaflet_frequency_low`: Vector of frequency values for the -1 leaflet types (low) along the rachis sub-sections..

Note that the length of the two vectors must be the same. It will define how many sub-sections the rachis is divided into
for this computation.

# Returns

A vector of NamedTuples representing the `(;high, medium, low)` frequencies for each sub-section.
"""
function compute_leaflet_type_frequencies(leaflet_frequency_high, leaflet_frequency_low)
    @assert length(leaflet_frequency_high) == length(leaflet_frequency_low) "Vectors must be of the same length"

    n = length(leaflet_frequency_high)
    leaflet_type_frequencies = Vector{NamedTuple{(:high, :medium, :low),Tuple{Float64,Float64,Float64}}}(undef, n)

    for i in 1:n
        medium_frequency = 1.0 - leaflet_frequency_high[i] - leaflet_frequency_low[i]
        @assert medium_frequency >= 0 "The sum of frequencies for high (+1) and low (-1) leaflets must be less than or equal to 1 for each section: section $i has a sum of $(leaflet_frequency_high[i] + leaflet_frequency_low[i])"
        leaflet_type_frequencies[i] = (high=leaflet_frequency_high[i], medium=medium_frequency, low=leaflet_frequency_low[i])
    end

    return leaflet_type_frequencies
end

"""
    draw_group_size(rel_pos, leaflet_type_frequencies, rng)

Determine the size of a leaflet group based on the relative position along the rachis and frequency patterns.

# Arguments

- `rel_pos`: Relative position on the rachis (0 to 1], where 0 is the base (C-point) and 1 is the tip (A-point).
- `leaflet_type_frequencies`: Vector of NamedTuples representing frequency distributions for each rachis segment, with fields:
  - `high`: Frequency of plane=+1 leaflets (first leaflet in each group), *i.e.* leaflets on "high" position
  - `medium`: Frequency of plane=0 leaflets (intermediate leaflets in groups), *i.e.* leaflets on "medium" position, horizontally inserted on the rachis
  - `low`: Frequency of plane=-1 leaflets (terminal leaflets in groups), *i.e.* leaflets on "low" position
- `rng`: Random number generator for stochastic determination.

# Details

This function implements an inverse relationship between the frequency of high (plane=1) leaflets 
and group size, modeling a fundamental biological pattern in palm frond architecture:

- Segments with high frequency of high leaflets produce many small groups of leaflets
- Segments with low frequency of high leaflets produce fewer, larger groups of leaflets

The calculation uses a probabilistic rounding mechanism to ensure proper statistical distribution 
of group sizes. This creates the natural variation in leaflet grouping patterns seen along real palm
fronds, where clustering patterns change systematically from base to tip.

# Returns

An integer representing the number of leaflets in the group.
"""
function draw_group_size(rel_pos, leaflet_type_frequencies, rng)
    segment = floor(Int, rel_pos * length(leaflet_type_frequencies))
    segment = clamp(segment, 1, length(leaflet_type_frequencies)) # Ensure segment is within bounds
    size_d = 1.0 / leaflet_type_frequencies[segment].high

    size_i = floor(Int, size_d)
    delta = size_d - size_i

    if rand(rng) < delta
        size_i += 1
    end

    return size_i
end