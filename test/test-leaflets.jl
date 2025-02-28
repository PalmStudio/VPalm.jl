
# Defining some parameter values;
parameters = Dict(
    "leaflets_nb_max" => 171.716496884477,
    "leaflets_nb_min" => 20,
    "leaflets_nb_slope" => 0.252868691666823,
    "leaflets_nb_inflexion" => 2.3256241328755,
    "nbLeaflets_SDP" => 0,
    "leaflet_position_shape_coefficient" => 2.47840369933858,
    "leaflet_frequency_high" => # Frequency of leaflets of type "+1" (High) along 10 sub-sections of the rachis
        [
            0.219712525667351,
            0.388663967611336,
            0.333333333333333,
            0.333333333333333,
            0.367647058823529,
            0.369158878504673,
            0.313084112149533,
            0.43801652892562,
            0.341176470588235,
            0.0688405797101449,
        ],
    "leaflet_frequency_low" => # Frequency of leaflets of type "-1" (Low) along 10 sub-sections of the rachis
        [
            0.0184804928131417,
            0.165991902834008,
            0.28310502283105,
            0.318627450980392,
            0.357843137254902,
            0.41588785046729,
            0.378504672897196,
            0.227272727272727,
            0.2,
            0.036231884057971,
        ],
    "rng" => Random.MersenneTwister(1234)
)

# Define The rachis length:
rachis_length = 2.0

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
    leaflets_grouped = VPalm.group_leaflets(leaflets_relative_positions, leaflets_type_frequency, parameters["rng"])
    @test leaflets_grouped == (group=[1, 1, 1, 1, 1, 2, 2, 2, 3, 3], group_size=[5, 5, 5, 5, 5, 3, 3, 3, 2, 2], plane=[1, -1, 0, 0, 0, 1, 0, 0, 1, 0]) == (group=[1, 1, 1, 1, 2, 2, 3, 3, 4, 4], group_size=[4, 4, 4, 4, 2, 2, 2, 2, 2, 2], plane=[1, 0, -1, 0, 1, -1, 1, 0, 1, 0])
end




leaflets = VPalm.group_leaflets(leaflets_relative_positions, leaflets_type_frequency, rng)
leaflets_relative_positions = VPalm.relative_leaflet_position.(collect(1:nb_leaflets) ./ nb_leaflets, parameters["leaflet_position_shape_coefficient"])
positions = leaflets_relative_positions .* rachis_length
VPalm.shrink_leaflets_in_groups!(positions, leaflets)
VPalm.normalize_positions!(initial_leaflets_position, rachis_length)




function test1(positions, leaflets)
    pos = deepcopy(positions)
    leaflets = deepcopy(leaflets)
    shrink_leaflets_in_groups!(pos, leaflets)

    return pos
end

function test2(positions, leaflets)
    pos = deepcopy(positions)
    leaflets = deepcopy(leaflets)
    apply_group_spacing!(pos, leaflets)

    return pos
end

@btime test1($positions, $leaflets);
@btime test2($positions, $leaflets);