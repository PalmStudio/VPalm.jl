"""
    interp_points(mat_points, pas)

Interpolate points along a curve.

# Arguments
- `mat_points`: Matrix of points (x, y, z coordinates).
- `pas`: Distance between interpolated points (m).

# Returns
- A NamedTuple with fields:
  - `points`: Meshes.Points of the interpolated points.
  - `i_discret_pts_exp`: Vector of indices of experimental points in the discretized curve.
  - `vec_dist_p2p1`: Vector of distances between consecutive interpolated points.
  - `vec_angle_xy`: Vector of angles between segments and the XY plane.
  - `vec_angle_xz`: Vector of angles between segments and the XZ plane.
"""
function interp_points(points, pas)
    # Distance and angles of each segment P2P1
    vdist_p2p1, vangle_xy, vangle_xz = xyz_to_dist_and_angles(points)

    dist_lineique = cumsum(vdist_p2p1)
    dist_totale = last(dist_lineique)

    # The distances of the segments cannot be zero
    # The origin point (0,0,0) cannot be in mat_points
    if any(vdist_p2p1 .== zero(vdist_p2p1[1]))
        error("Found distances between segments equal to 0.")
    end

    # Construction of interpolated points
    nlin = round(Int, dist_totale / pas) + 1
    pas = dist_totale / (nlin - 1)

    npoints_exp = length(points)

    vec_dist = [zero(eltype(pas)); fill(pas, nlin - 1)]
    dist_interp = cumsum(vec_dist)

    # Avoid rounding errors
    dist_lineique[npoints_exp] = dist_lineique[npoints_exp] + (1.0 * unit(eltype(dist_lineique)))

    mat_xyz = zeros(3, 0)

    for i in 1:npoints_exp
        if i == 1
            ind_points = (dist_interp .<= dist_lineique[i])
        else
            ind_points = (dist_interp .> dist_lineique[i-1]) .& (dist_interp .<= dist_lineique[i])
        end

        if !any(ind_points)
            error("No point found")
        end

        dist_points = vec_dist[ind_points]
        vec_rot = RotYZ(-vangle_xy[i], -vangle_xz[i]) * [dist_points'; zeros(eltype(dist_points), 2, length(dist_points))]

        if i > 1
            vec_points = cumsum(vec_rot, dims=2) .+ mat_xyz[:, end]
            mat_xyz = hcat(mat_xyz, vec_points)
        else
            mat_xyz = cumsum(vec_rot, dims=2)
        end
    end

    vec_x = mat_xyz[1, :]
    vec_y = mat_xyz[2, :]
    vec_z = mat_xyz[3, :]

    # Identification of experimental points
    # in the linear discretization
    i_discret_pts_exp = zeros(Int, npoints_exp)

    for i in 1:npoints_exp
        p = Meshes.coords(points[i])
        equad = sqrt.((vec_x .- p.x) .^ 2 .+ (vec_y .- p.y) .^ 2 .+ (vec_z .- p.z) .^ 2)
        ind = findall(equad .== minimum(equad))
        i_discret_pts_exp[i] = ind[1]
    end

    # Make points out the coordinates:
    interp_points = [Meshes.Point(vec_x[i], vec_y[i], vec_z[i]) for i in eachindex(vec_x)]

    # Distance and angles of the interpolated points
    XYZ_Agl = xyz_to_dist_and_angles(interp_points)

    return (points=interp_points, i_discret_pts_exp=i_discret_pts_exp, vec_dist_p2p1=XYZ_Agl.dist_p2p1,
        vec_angle_xy=XYZ_Agl.vangle_xy, vec_angle_xz=XYZ_Agl.vangle_xz)
end