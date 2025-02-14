# Define dummy functions if they are not yet implemented
# Create a dummy test file
df = CSV.read(joinpath(@__DIR__, "files/6_EW01.22_17_kanan_unbent.csv"), DataFrame)

# Set test parameters
pas = 0.02  # in meter. -> Length of the segments that discretize the object.
Ncalc = 100 # number of points used in the grid that discretized the section.
Nboucle = 15 # if we want to compute the torsion after the bending step by step instead of
elastic_modulus = 2000.0
shear_modulus = 400.0

ref = CSV.read(joinpath(@__DIR__, "files/6_EW01.22_17_kanan_unbent_bend.csv"), DataFrame)
@testset "bend works" begin
    # Dummy input data for bend function
    # Call the function
    out = VPalm.bend(df.type, df.width, df.height, df.torsion, df.x, df.y, df.z, df.mass, df.mass_right, df.mass_left,
        df.distance_application, elastic_modulus, shear_modulus, pas, Ncalc, Nboucle)
    # CSV.write(joinpath(@__DIR__, "files/6_EW01.22_17_kanan_unbent_bend.csv"), DataFrame(out))
    @test isapprox(ref, DataFrame(out); atol=10)
end

