using VPalm
using Test
using Aqua
using JET

@testset "VPalm.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(VPalm)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(VPalm; target_defined_modules = true)
    end
    # Write your tests here.
end
