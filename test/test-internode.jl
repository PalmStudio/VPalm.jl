file = joinpath(dirname(dirname(pathof(VPalm))), "test", "files", "parameter_file.yml")
parameters = read_parameters(file)

@testset "internode diameter" begin
end

@testset "internode length" begin
end

@testset "phyllotactic angle" begin
    @test VPalm.phyllotactic_angle(0.0, 0.05, 0) == 0.0
    @test VPalm.phyllotactic_angle(45.0, 0.05, 0) == 45.0
    @test VPalm.phyllotactic_angle(45.0, 1.0) != 45.0
end