file = joinpath(dirname(dirname(pathof(VPalm))), "test", "files", "parameter_file.yml")
parameters = read_parameters(file)

@testset "static mockup" begin
    mtg = VPalm.mtg_skeleton(parameters)
    @test length(mtg) == (parameters["nb_leaves_emitted"] + parameters["nb_internodes_before_planting"]) * 3 + 2 # nb leaves emitted * 3 (phytomer + internode + leaf) + 2 (stem + plant)
    # Check the length of the mockup: nb leaves emitted * 3 (phytomer + internode + leaf) + 2 (stem + plant)
    @test typeof(mtg) == MultiScaleTreeGraph.Node{MultiScaleTreeGraph.NodeMTG,Dict{Symbol,Any}}
    @test mtg[1][:stem_bending] == 0.0
end
