@testset "dist_and_angles_to_xyz" begin
    # Test different vector sizes
    @test_throws ErrorException VPalm.dist_and_angles_to_xyz([1.0], [1.0, 2.0], [1.0])
    @test_throws ErrorException VPalm.dist_and_angles_to_xyz([1.0], [1.0], [1.0, 2.0])

    # Test simple case (90° to XY, 0° to XZ)
    dist = [1.0]
    angle_xy = [π / 2] # 90°
    angle_xz = [0.0]
    res = VPalm.dist_and_angles_to_xyz([1.0], angle_xy, angle_xz)
    @test isapprox(res.vec_x, [0.0]; atol=1e-10)
    @test isapprox(res.vec_y, [0.0]; atol=1e-10)
    @test isapprox(res.vec_z, [1.0]; atol=1e-10)

    # Test simple case (0° to XY, 0° to XZ)
    dist = [1.0]
    angle_xy = [0.0]
    angle_xz = [0.0]
    res = VPalm.dist_and_angles_to_xyz(dist, angle_xy, angle_xz)
    @test isapprox(res.vec_x, [1.0]; atol=1e-10)
    @test isapprox(res.vec_y, [0.0]; atol=1e-10)
    @test isapprox(res.vec_z, [0.0]; atol=1e-10)

    # Test 2 segments
    dist = [1.0, 1.0]
    angle_xy = [0.0, π / 4]  # 0° then 45°
    angle_xz = [0.0, 0.0]
    res = VPalm.dist_and_angles_to_xyz(dist, angle_xy, angle_xz)

    # First point
    @test isapprox(res.vec_x[1], 1.0; atol=1e-10)
    @test isapprox(res.vec_y[1], 0.0; atol=1e-10)
    @test isapprox(res.vec_z[1], 0.0; atol=1e-10)

    # Second point
    @test isapprox(res.vec_x[2], 1.0 + cos(π / 4); atol=1e-10)
    @test isapprox(res.vec_y[2], 0.0; atol=1e-10)
    @test isapprox(res.vec_z[2], sin(π / 4); atol=1e-10)

    # Test with 90° in the XZ plane rotation
    dist = [1.0]
    angles_xy = [0.0]
    angles_xz = [π / 2]  # 90° dans le plan XZ
    res = VPalm.dist_and_angles_to_xyz(dist, angles_xy, angles_xz)
    @test isapprox(res.vec_x[1]; 0.0; atol=1e-10)
    @test isapprox(res.vec_y[1]; 1.0; atol=1e-10)
    @test isapprox(res.vec_z[1]; 0.0; atol=1e-10)
end


@testitem "xyz_to_dist_and_angles" begin
    # Test different vector sizes
    @test_throws ErrorException VPalm.xyz_to_dist_and_angles([1.0], [1.0, 2.0], [1.0])
    @test_throws ErrorException VPalm.xyz_to_dist_and_angles([1.0], [1.0], [1.0, 2.0])

    # Test simple case (90° to XY, 0° to XZ)
    x = [0.0]
    y = [0.0]
    z = [1.0]
    res = VPalm.xyz_to_dist_and_angles(x, y, z)
    @test isapprox(res.dist_p2p1, [1.0]; atol=1e-10)
    @test isapprox(res.vangle_xy, [π / 2]; atol=1e-10)
    @test isapprox(res.vangle_xz, [0.0]; atol=1e-10)

    # Test simple case (0° to XY, 0° to XZ)
    x = [1.0]
    y = [0.0]
    z = [0.0]
    res = VPalm.xyz_to_dist_and_angles(x, y, z)
    @test isapprox(res.dist_p2p1, [1.0]; atol=1e-10)
    @test isapprox(res.vangle_xy, [0.0]; atol=1e-10)
    @test isapprox(res.vangle_xz, [0.0]; atol=1e-10)

    # Test 2 segments
    x = [1.0, 1.0 + cos(π / 4)]
    y = [0.0, 0.0]
    z = [0.0, sin(π / 4)]
    res = VPalm.xyz_to_dist_and_angles(x, y, z)

    # First point
    @test isapprox(res.dist_p2p1[1], 1.0; atol=1e-10)
    @test isapprox(res.vangle_xy[1], 0.0; atol=1e-10)
    @test isapprox(res.vangle_xz[1], 0.0; atol=1e-10)

    # Second point
    @test isapprox(res.dist_p2p1[2], 1
        ; atol=1e-10)
    @test isapprox(res.vangle_xy[2], π / 4; atol=1e-10)
    @test isapprox(res.vangle_xz[2], 0.0; atol=1e-10)
end