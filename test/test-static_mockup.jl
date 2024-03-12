file = joinpath(dirname(dirname(pathof(VPalm))), "test", "files", "parameter_file.yml")

@testset "static mockup" begin
    parameters = read_parameters(file)

    mockup = static_mockup(parameters)
    # Check the length of the mockup: nb leaves emitted * 3 (phytomer + internode + leaf) + 2 (stem + plant)
    @test typeof(mockup) == MultiScaleTreeGraph.Node{MultiScaleTreeGraph.NodeMTG,Dict{Symbol,Any}}
    @test length(mockup) == (parameters["nb_leaves_emitted"]) * 3 + 2
end
