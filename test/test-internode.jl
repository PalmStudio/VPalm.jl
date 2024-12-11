file = joinpath(dirname(dirname(pathof(VPalm))), "test", "files", "parameter_file.yml")
parameters = read_parameters(file)

@testset "internode diameter" begin
end

@testset "internode length" begin
    @test VPalm.internode_length(1, 120, 3, 9, 20, 0.001) ≈ 0.001
    @test VPalm.internode_length(20, 120, 3, 9, 20, 0.001) ≈ 0.02388
    @test VPalm.internode_length(21, 120, 3, 9, 20, 0.001) ≈ 0.02388
    @test VPalm.internode_length(131, 120, 3, 9, 20, 0.001) ≈ 0.02388
    @test VPalm.internode_length(140, 120, 3, 9, 20, 0.001) ≈ 0.001
    @test sum(VPalm.internode_length.(1:140, 120, 3, 9, 20, 0.001)) ≈ 3
end

@testset "phyllotactic angle" begin
    @test VPalm.phyllotactic_angle(0.0, 0.05, 0) == 0.0
    @test VPalm.phyllotactic_angle(45.0, 0.05, 0) == 45.0
    @test VPalm.phyllotactic_angle(45.0, 1.0) != 45.0
end