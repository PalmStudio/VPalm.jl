
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

    create_single_leaflet(
        Ref(1),
        1,              # Index 1
        1,              # Scale 2 (typical for leaflets)
        10,
        0.5,
        0.5,
        1,
        1,
        1.0,
        0.2,
        offset=0.0,  # No offset needed for standalone
        rng=rng
    )
end