file = joinpath(dirname(dirname(pathof(VPalm))), "test", "files", "parameter_file.yml")

@testset "read_parameters" begin
    params = read_parameters(file)
    @test params["seed"] == 0
    @test params["frondTipHeight"] == 0.5
    @test params["rachisBiomass"] == [2607.60521189879, 2582.76405648725, 2557.92290107571, 2533.08174566417, 2508.24059025263]
end

@testset "write_parameters" begin
    params = read_parameters(file)

    params2 = mktemp() do f, io
        write_parameters(f, params)
        params2 = read_parameters(f)
        return params2
    end

    @test params == params2
end