file = joinpath(dirname(dirname(pathof(VPalm))), "test", "files", "parameter_file.yml")
parameters = read_parameters(file)

@testset "stem bending" begin
    @test VPalm.stem_bending(parameters["stem_bending_mean"], parameters["stem_bending_sd"]) == 0.0
    @test VPalm.stem_bending(45.0, parameters["stem_bending_sd"]) == 45.0
    @test VPalm.stem_bending(45.0, 1.0) != 45.0
end

@testset "stem height" begin
    # Before stem_growth_start:
    @test VPalm.stem_height(100, parameters["initial_stem_height"], parameters["stem_height_coefficient"], parameters["internode_length_at_maturity"], 120.0, 0.0) == 0.3005242665305668

    # After stem_growth_start:
    @test VPalm.stem_height(parameters["nb_leaves_emitted"], parameters["initial_stem_height"], parameters["stem_height_coefficient"], parameters["internode_length_at_maturity"], parameters["stem_growth_start"], 0.0) == 1.1801911326167471
end

@testset "stem diameter" begin
    @test VPalm.stem_diameter(5.0335230597,
    parameters["stem_diameter_max"],
    parameters["stem_diameter_slope"],
    parameters["stem_diameter_inflection"],
    parameters["stem_diameter_residual"],
    parameters["stem_diameter_snag"],
    0.0) == 0.16228677550579518