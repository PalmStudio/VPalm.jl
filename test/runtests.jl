using VPalm
using Test
using Aqua
using JET

@testset "VPalm.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(VPalm, ambiguities=false)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(VPalm; target_defined_modules=true)
    end

    @testset "Parameters IO" begin
        include("test-parameters_IO.jl")
    end

    @testset "Static mockup" begin
        include("test-static_mockup.jl")
    end
    # Write your tests here.
end
