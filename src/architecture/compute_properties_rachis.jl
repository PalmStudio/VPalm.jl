"""
    biomechanical_properties_rachis(
        rachis_twist_initial_angle, rachis_twist_initial_angle_sdp,
        elastic_modulus, shear_modulus, rachis_length,
        lenflet_length_at_b_intercept, leaflet_length_at_b_slope, relative_position_bpoint,
        relative_position_bpoint_sd, relative_length_first_leaflet, relative_length_last_leaflet, relative_position_leaflet_max_length,
        rachis_fresh_weigth, rank, height_cpoint, zenithal_cpoint_angle, nb_sections,
        height_rachis_tappering,
        points, iterations, angle_max,
        rng
    )

Use of the biomechanical model to compute the properties of the rachis.

# Arguments
- `rachis_twist_initial_angle`: initial twist angle of the rachis (°)
- `rachis_twist_initial_angle_sdp`: standard deviation of the initial twist angle of the rachis (°)
- `elastic_modulus`: elastic modulus of the rachis (Pa)
- `shear_modulus`: shear modulus of the rachis (Pa)
- `rachis_length`: length of the rachis (m)
- `lenflet_length_at_b_intercept`: intercept of the linear function for the leaflet length at the B point (m)
- `leaflet_length_at_b_slope`: slope of the linear function for the leaflet length at the B point (m)
- `relative_position_bpoint`: relative position of the B point on the rachis (0: base to 1: tip)
- `relative_position_bpoint_sd`: standard deviation of the relative position of the B point on the rachis
- `relative_length_first_leaflet`: relative length of the first leaflet on the rachis (0 to 1)
- `relative_length_last_leaflet`: relative length of the last leaflet on the rachis (0 to 1)
- `relative_position_leaflet_max_length`: relative position of the longest leaflet on the rachis (0.111 to 0.999)
- `rachis_fresh_weigth`: fresh weight of the rachis (g)
- `rank`: rank of the rachis
- `height_cpoint`: height of the C point (m)
- `zenithal_cpoint_angle`: zenithal angle of the C point (°)
- `nb_sections`: number of sections to compute the bending
- `height_rachis_tappering`: tappering factor for the rachis height
- `points`: number of points to compute the bending
- `iterations`: number of iterations to compute the bending
- `angle_max`: maximum angle to compute the bending (°)
- `rng`: the random number generator

# Returns
A named tuple with the following fields:
- `length`: vector with the length of each segment
- `points_positions`: the position of the points along the rachis
- `bending`: the bending angle of the rachis
- `deviation`: the deviation of the rachis (angle in the xz plane)
- `torsion`: the torsion of the rachis
- `x`: the x coordinates of the rachis
- `y`: the y coordinates of the rachis
- `z`: the z coordinates of the rachis

# Details
Split the rachis into 5 segments defined by remarkable points (C, C-B, B, B-A, A).
Each segment has a particular shape, a mass, and the leaflets on both sides of the rachis have a mass.
Coefficents are used to compute the mass distribution and relative lengths of segments.
The rachis is bent using the `bend` function.
"""
function biomechanical_properties_rachis(
    rachis_twist_initial_angle, rachis_twist_initial_angle_sdp,
    elastic_modulus, shear_modulus, rachis_length,
    lenflet_length_at_b_intercept, leaflet_length_at_b_slope, relative_position_bpoint,
    relative_position_bpoint_sd, relative_length_first_leaflet, relative_length_last_leaflet, relative_position_leaflet_max_length,
    rachis_fresh_weigth, rank, height_cpoint, zenithal_cpoint_angle, nb_sections,
    height_rachis_tappering,
    points, iterations, angle_max,
    rng
)
    # Frond section types (e.g., rectangle, ellipsoid, etc.)
    type = [1, 2, 3, 4, 5]
    npoints = length(type)

    # Compute initial torsion (using prms and rnd assumed to be defined)
    initial_torsion_sdp = rachis_twist_initial_angle + normal_deviation_draw(rachis_twist_initial_angle_sdp, rng)

    initial_torsion_vec = fill(initial_torsion_sdp, npoints)
    # Relative position of the remarkable points (C, C-B, B, B-A, A) on the rachis:
    relative_position_remarkable_points = [0.0000001, 0.336231351023383, 0.672462702046766, 0.836231351023383, 1.0]
    # Note: we use those positions as remarkable points along the rachis, and each segment (or section) is defined by two consecutive points.
    # Each segment has a particular shape, a mass, and the leaflets on both sides of the rachis have a mass.

    # Relative position at the middle of each segment:
    relative_position_mid_segment = [0.1681157, 0.504347, 0.672462702046766, 0.754347, 0.9181157]

    # Distribution of the mass for each segment relative to the total rachis mass:
    mass_distribution_segment_rachis = [0.0, 0.648524097435024, 0.277401814695433, 0.0601164171693578, 0.0139576707001849]

    # Distribution of the mass for each leaflet relative to the total rachis mass:
    mass_distribution_segment_leaflet = [0.0, 0.0658151279405379, 0.201957451540734, 0.105263443497354, 0.0475385258600695]

    # Initialization of data computed for each of the 5 remarkable points:
    mass = zeros(Float64, npoints)               # Mass of each segment represented by the points
    mass_right = zeros(Float64, npoints)         # Mass of the leaflets on the right-hand side of each segment
    mass_left = zeros(Float64, npoints)          # Mass of the leaflets on the left-hand side of each segment
    width_bend = zeros(Float64, npoints)         # Width of the segment (rachis width)
    height_bend = zeros(Float64, npoints)        # Height of the segment (rachis height)
    distances = zeros(Float64, npoints)          # Distance between the points projected on the X axis
    distance_application = zeros(Float64, npoints) # Application distance for forces (if needed)

    leaflet_length_at_bpoint = length_at_bpoint(rachis_length, lenflet_length_at_b_intercept, leaflet_length_at_b_slope)
    leafletLengthMax = max_leaflet_length(leaflet_length_at_bpoint, relative_position_bpoint, relative_position_bpoint_sd, relative_length_first_leaflet, relative_length_last_leaflet, relative_position_leaflet_max_length, rng)

    # Parameters to compute rachis width from rachis height:
    ratioPointC = 0.5220
    ratioPointA = 1.0053
    posRatioMax = 0.6636
    ratioMax = 1.5789

    for i in 1:npoints
        distances[i] = rachis_length * relative_position_remarkable_points[i]
        mass[i] = mass_distribution_segment_rachis[i] * rachis_fresh_weigth
        # we consider that the leaflets on both sides of the rachis have the same mass:
        mass_right[i] = mass_distribution_segment_leaflet[i] * rachis_fresh_weigth
        mass_left[i] = mass_distribution_segment_leaflet[i] * rachis_fresh_weigth

        # leaflet length at the middle of the segment (in m):
        length_leaflets_segment = leafletLengthMax * relative_leaflet_length(
            relative_position_mid_segment[i],
            relative_length_first_leaflet, relative_length_last_leaflet,
            relative_position_leaflet_max_length
        )

        distance_application[i] = length_leaflets_segment / 10.0  # The leaflet weight is applied at the middle of the leaflets

        if rank < 3
            distance_application[i] = 1e-8
            initial_torsion_vec[i] = 0.0
        end

        height_bend[i] = rachis_height(relative_position_remarkable_points[i], height_cpoint, height_rachis_tappering)
        width_bend[i] = height_bend[i] / height_to_width_ratio(relative_position_remarkable_points[i], ratioPointC, ratioPointA, posRatioMax, ratioMax)
    end

    # Un-bent coordinates (take the leaf as a straight line in x and z)
    x = zeros(5)
    y = zeros(5)
    z = zeros(5)

    for n in eachindex(distances)
        x[n] = cosd(zenithal_cpoint_angle) * distances[n] # Note: zenithal_cpoint_angle is in degrees, so we use cosd instead of cos
        z[n] = sind(zenithal_cpoint_angle) * distances[n]
    end

    step = rachis_length / (nb_sections - 1)
    # @show type width_bend height_bend initial_torsion_vec x y z mass mass_right mass_left distance_application elastic_modulus shear_modulus
    # error("stop here")

    # Call the bend function, which returns a vector of arrays:
    # bending -> { PtsX, PtsY, PtsZ, PtsDist, PtsAglXY, PtsAglXZ, PtsAglTor }
    bending = bend(
        type, width_bend, height_bend, initial_torsion_vec, x, y, z, mass, mass_right, mass_left,
        distance_application, elastic_modulus, shear_modulus, step, points, iterations;
        verbose=false, all_points=true, angle_max=angle_max
    )

    points_bending = .-bending.angle_xy
    points_bending[1] = -zenithal_cpoint_angle         # Initialize the first angle as the angle at C point

    return (
        length=fill(step, length(bending.x)), points_positions=bending.length, bending=points_bending, deviation=bending.angle_xz, torsion=bending.torsion,
        x=bending.x, y=bending.y, z=bending.z
    )
end


function length_at_bpoint(rachis_length, lenflet_length_at_b_intercept, leaflet_length_at_b_slope)
    return linear(rachis_length, lenflet_length_at_b_intercept, leaflet_length_at_b_slope)
end

function max_leaflet_length(leaflet_length_at_bpoint, relative_position_bpoint, relative_position_bpoint_sd, relative_length_first_leaflet, relative_length_last_leaflet, relative_position_leaflet_max_length, rng)
    relative_position_bpoint = relative_position_bpoint + normal_deviation_draw(relative_position_bpoint_sd, rng)
    return leaflet_length_at_bpoint / relative_leaflet_length(relative_position_bpoint, relative_length_first_leaflet, relative_length_last_leaflet, relative_position_leaflet_max_length)
end


"""
    relative_leaflet_length(x, relative_length_first_leaflet, relative_length_last_leaflet, relative_position_leaflet_max_length)

Relative leaflet length given by their relative position along the rachis.

# Arguments

- `x`: relative leaflet position on the rachis (0: base to 1: tip)
- `relative_length_first_leaflet`: relative length of the first leaflet on the rachis (0 to 1)
- `relative_length_last_leaflet`: relative length of the last leaflet on the rachis  (0 to 1)
- `relative_position_leaflet_max_length`: relative position of the longest leaflet on the rachis (0.111 to 0.999)
"""
function relative_leaflet_length(x, relative_length_first_leaflet, relative_length_last_leaflet, relative_position_leaflet_max_length)
    if x < relative_position_leaflet_max_length
        return relative_length_first_leaflet + ((1 - relative_length_first_leaflet) * x * (2 * relative_position_leaflet_max_length - x)) / relative_position_leaflet_max_length^2
    else
        return 1 + (relative_length_last_leaflet - 1) * (x - relative_position_leaflet_max_length)^2 / (1 - relative_position_leaflet_max_length)^2
    end
end


"""
    rachis_height(relative_position, cpoint_height, rachis_height_tappering)

Computes the rachis height (m) at a given relative position using a the height at C Point and rachis tappering.

# Arguments

- `relative_position`: The relative position along the rachis (0: base to 1: tip).
- `cpoint_height`: The height of the rachis at the C point, *i.e.* rachis base (m).
- `rachis_height_tappering`: The tappering factor for the rachis height.
"""
function rachis_height(relative_position, cpoint_height, rachis_height_tappering)
    return (1.0 + rachis_height_tappering * (relative_position^3)) * cpoint_height
end



"""
    rachis_width(relative_position, cpoint_width, rachis_width_tip)

Computes the rachis width (m) at a given relative position using the width at C Point and rachis width at the tip.

# Arguments

- `relative_position`: The relative position along the rachis (0: base to 1: tip).
- `cpoint_width`: The width of the rachis at the C point, *i.e.* rachis base (m).
- `rachis_width_tip`: The width of the rachis at the tip (m).
"""
function rachis_width(relative_position, cpoint_width, rachis_width_tip)
    return cpoint_width * (1.0 - relative_position) + rachis_width_tip * relative_position
end

"""
    height_to_width_ratio(x, ratio_point_c, ratio_point_a, pos_ratio_max, ratio_max)

Computes the relative width along the rachis.

# Arguments

- `x`: relative position on the rachis
- `ratio_point_c`: ratio at point C
- `ratio_point_a`: ratio at point A
- `pos_ratio_max`: relative position of the maximum value of the ratio
- `ratio_max`: maximum ratio value
"""
function height_to_width_ratio(x, ratio_point_c, ratio_point_a, pos_ratio_max, ratio_max)
    if x < pos_ratio_max
        return ratio_point_c + x * (ratio_max - ratio_point_c) / pos_ratio_max
    else
        return x * (ratio_point_a - ratio_max) / (1 - pos_ratio_max) +
               (ratio_max - pos_ratio_max * ratio_point_a) / (1 - pos_ratio_max)
    end
end