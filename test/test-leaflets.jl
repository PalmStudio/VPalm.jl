
# Defining some parameter values;
parameters = read_parameters(joinpath(dirname(dirname(pathof(VPalm))), "test", "files", "parameter_file.yml"))
# Define The rachis length:
rachis_length = 2.0u"m"
rng = Random.MersenneTwister(parameters["seed"])
nb_leaflets = VPalm.compute_number_of_leaflets(rachis_length, parameters["leaflets_nb_max"], parameters["leaflets_nb_min"], parameters["leaflets_nb_slope"], parameters["leaflets_nb_inflexion"], parameters["nbLeaflets_SDP"]; rng=rng)

@testset "Number of leaflets" begin
    @test nb_leaflets == 72
end

leaflets_type_frequency = VPalm.compute_leaflet_type_frequencies(parameters["leaflet_frequency_high"], parameters["leaflet_frequency_low"])
@testset "Leaflets type frequency" begin
    @test length(leaflets_type_frequency) == 10
    @test map(x -> x.high, leaflets_type_frequency) == parameters["leaflet_frequency_high"]
    @test map(x -> x.low, leaflets_type_frequency) == parameters["leaflet_frequency_low"]
    @test leaflets_type_frequency[1].medium == 0.7618069815195073
    @test leaflets_type_frequency[end].medium == 0.8949275362318841
end

leaflets_relative_positions = VPalm.relative_leaflet_position.(collect(1:nb_leaflets) ./ nb_leaflets, parameters["leaflet_position_shape_coefficient"])
@testset "Relative leaflet positions" begin
    @test length(leaflets_relative_positions) == nb_leaflets
    @test leaflets_relative_positions[1] == 0.0006706677303153433
    @test leaflets_relative_positions[end] == 1.0
    @test leaflets_relative_positions[35] == 0.5183724388916218
end


# Testing grouping on 10 leaflets only to make tests more readable:
@testset "Grouping leaflets" begin
    leaflets_relative_positions = VPalm.relative_leaflet_position.(collect(1:10) ./ 10, parameters["leaflet_position_shape_coefficient"])
    leaflets_grouped = VPalm.group_leaflets(leaflets_relative_positions, leaflets_type_frequency, rng)
    groups = []
    i = Ref(1)
    group = Ref(1)
    while true
        if i[] > length(leaflets_grouped.group_size)
            break
        end
        append!(groups, fill(group[], leaflets_grouped.group_size[i[]]))
        i[] += leaflets_grouped.group_size[i[]]
        group[] += 1
    end
    @test leaflets_grouped.group == groups
end

@testset "Leaflets dimensions" begin
    leaflet_length_at_b = VPalm.leaflet_length_at_bpoint(rachis_length, parameters["leaflet_length_at_b_intercept"], parameters["leaflet_length_at_b_slope"])
    @test leaflet_length_at_b ≈ 0.7192367735814625u"m"
    leaflet_max_length = VPalm.leaflet_length_max(leaflet_length_at_b, parameters["relative_position_bpoint"], parameters["relative_length_first_leaflet"], parameters["relative_length_last_leaflet"], parameters["relative_position_leaflet_max_length"])
    @test leaflet_max_length ≈ 0.7483089246613868u"m"
    leaflet_width_at_b = VPalm.leaflet_width_at_bpoint(rachis_length, parameters["leaflet_width_at_b_intercept"], parameters["leaflet_width_at_b_slope"])
    @test leaflet_width_at_b ≈ 0.05501868507830992u"m"
    leaflet_max_width = VPalm.leaflet_width_max(leaflet_width_at_b, parameters["relative_position_bpoint"], parameters["relative_width_first_leaflet"], parameters["relative_width_last_leaflet"], parameters["relative_position_leaflet_max_width"])
    @test leaflet_max_width ≈ 0.05821782360344443u"m"
end

@testset "Make a leaflet" begin
    parameters_ = copy(parameters)
    parameters_["leaflet_stiffness_sd"] = 0.0
    parameters_["leaflet_axial_angle_sdp"] = 0.0
    leaflet = VPalm.create_single_leaflet(
        Ref(1),
        1,              # Index 1
        1,              # Scale 2 (typical for leaflets)
        10,             # Rank 10
        0.5,            # Relative position
        0.5,            # Normalized leaflet rank
        1,              # Plane (low:-1, planar: 0 or high:1)
        1,              # Side of the rachis (1 or -1)
        0.7483089246613868u"m", # Leaflet maximum length
        0.05821782360344443u"m", # Leaflet maximum width
        parameters_,     # Parameters
        offset=0.0,     # No offset needed for standalone
        rng=rng
    )


    segment_position = leaflet[1].segment_boundaries * leaflet.length
    # Compute the length of each segment from each end position and previous segment end position:
    segment_length = diff([0.0u"m"; segment_position])
    @test sum(segment_length) ≈ leaflet.length
    @test sum(descendants(leaflet, :length)) ≈ leaflet.length

    @test length(leaflet) == 6 # One node for the leaflet, 5 for the leaflet segments
    @test symbol(leaflet) == "Leaflet"
    @test symbol(leaflet[1]) == "LeafletSegment"
    @test leaflet.leaflet_rank == 0.5
    @test leaflet.length ≈ 0.7471756197471424u"m"
    @test leaflet.offset == 0.0
    @test leaflet.plane == 1
    @test leaflet.relative_position == 0.5
    # @test leaflet.rot_bearer_x ≈ 38.10350146807654 # This one has randomness baked into it
    @test leaflet.rot_bearer_z ≈ 71.35581494643947
    @test leaflet.side == 1
    @test leaflet.width ≈ 0.04985922542234243u"m"
    @test leaflet.twist_angle == 10.0
    @test leaflet.stiffness ≈ 1500.0
    @test leaflet.tapering == 0.5
end