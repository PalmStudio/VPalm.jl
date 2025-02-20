using VPalm
using Test
using Aqua
using JET
using MultiScaleTreeGraph

@testset "VPalm.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(VPalm, ambiguities=false)
    end

    if VERSION >= v"1.10"
        # See this issue: https://github.com/aviatesk/JET.jl/issues/665
        @testset "Code linting (JET.jl)" begin
            JET.test_package(VPalm; target_defined_modules=true)
        end
    end

    @testset "Parameters IO" begin
        include("test-parameters_IO.jl")
    end

    @testset "Stem allometries" begin
        include("test-stem.jl")
    end

    @testset "Static mockup" begin
        include("test-static_mockup.jl")
    end
    # Write your tests here.
end
