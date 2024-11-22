file = joinpath(dirname(dirname(pathof(VPalm))), "test", "files", "6_EW01.22_17_kanan.txt")
df = read_mat(filepath)
dfu = unbend(df)

pas= 0.02 # in meter. -> Length of the segments that discretize the object.
Ncalc= 100 # number of points used in the grid that discretized the section.
Nboucle= 15 # if we want to compute the torsion after the bending step by step instead of

@testset "bend works" begin
    result = bend(data = dfu, elastic_modulus =  2000, shear_modulus = 400,
        pas, Ncalc, Nboucle, verbose = TRUE)
    expected_result = CSV.read("bend.test") # Assuming "bend.test" is a file with the expected result in a similar format
    @test result == expected_result
end

@testset "unbend works" begin
    @test round.(dfu.x, digits=5) == [0.00066, 0.8833, 1.76595, 2.27314, 2.78033]
    @test round.(dfu.y, digits=5) == [0, 0, 0, 0, 0]
    @test round.(dfu.z, digits=5) == [0.00075, 1.00899, 2.01722, 2.59658, 3.17594]
end