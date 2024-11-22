file = joinpath(dirname(dirname(pathof(VPalm))), "test", "files", "parameter_file.yml")
parameters = read_parameters(file)

@testset "static mockup" begin
    mtg = VPalm.mtg_skeleton(parameters["nb_leaves_emitted"])
    mtg[1]
    mockup = static_mockup(parameters)
    # Check the length of the mockup: nb leaves emitted * 3 (phytomer + internode + leaf) + 2 (stem + plant)
    @test typeof(mockup) == MultiScaleTreeGraph.Node{MultiScaleTreeGraph.NodeMTG,Dict{Symbol,Any}}
    @test length(mockup) == (parameters["nb_leaves_emitted"]) * 3 + 2
    @test mockup[1][:stem_bending] == 0.0
end
