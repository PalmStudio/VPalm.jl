"""
    xyz_to_dist_and_angles(vec_x, vec_y, vec_z)

Compute segment lengths and angles from point coordinates.

# Arguments

- `vec_x`: Vector of X coordinates (m).
- `vec_y`: Vector of Y coordinates (m).
- `vec_z`: Vector of Z coordinates (m).

# Returns

- A NamedTuple with fields:
  - `dist_p2p1`: Vector of segment lengths (m).
  - `vangle_xy`: Vector of angles between the segment and the XY plane (radians).
  - `vangle_xz`: Vector of angles between the segment and the XZ plane (radians).
"""
function xyz_to_dist_and_angles(vec_x, vec_y, vec_z)
    n = length(vec_x)

    if length(vec_y) != n
        error("Length of Y coordinates not equal to X coordinates")
    end
    if length(vec_z) != n
        error("Length of Z coordinates not equal to X coordinates")
    end

    dist_p2p1 = zeros(n)
    vangle_xy = zeros(n)
    vangle_xz = zeros(n)

    for iter in 1:n
        p2 = [vec_x[iter], vec_y[iter], vec_z[iter]]

        if iter == 1
            p1 = zeros(3)
        else
            p1 = [vec_x[iter-1], vec_y[iter-1], vec_z[iter-1]]
        end

        p2p1 = p2 .- p1

        if p2p1[1] != 0 #! why?? We could have some distance in the Z or Y axis?
            # Distances
            dist_p2p1[iter] = sqrt(p2p1[1]^2 + p2p1[2]^2 + p2p1[3]^2)

            # Angles
            vangle_xy[iter] = atan(p2p1[3] / sqrt(p2p1[1]^2 + p2p1[2]^2))
            vangle_xz[iter] = atan(p2p1[2] / p2p1[1])
        end
    end

    return (dist_p2p1=dist_p2p1, vangle_xy=vangle_xy, vangle_xz=vangle_xz)
end

"""
    dist_and_angles_to_xyz(dist_p2p1, vangle_xy, vangle_xz)

Transform distances and angles into point coordinates.

# Arguments
- `dist_p2p1`: Vector of segment lengths (m).
- `vangle_xy`: Vector of angles between the segment and the XY plane (radians).
- `vangle_xz`: Vector of angles between the segment and the XZ plane (radians).

# Returns
- A NamedTuple with fields:
  - `vec_x`: Vector of x coordinates.
  - `vec_y`: Vector of y coordinates.
  - `vec_z`: Vector of z coordinates.
"""
function dist_and_angles_to_xyz(dist_p2p1, vangle_xy, vangle_xz)
    n = length(dist_p2p1)

    if length(vangle_xy) != n
        error("length of vangle_xy != n")
    end
    if length(vangle_xz) != n
        error("length of vangle_xz != n")
    end

    vec_x = zeros(n)
    vec_y = zeros(n)
    vec_z = zeros(n)

    for iter in 1:n
        dz = dist_p2p1[iter] * sin(vangle_xy[iter])
        dist_xy = dist_p2p1[iter] * cos(vangle_xy[iter])

        dx = dist_xy * cos(vangle_xz[iter])
        dy = dist_xy * sin(vangle_xz[iter])

        if iter == 1
            vec_x[iter] = dx
            vec_y[iter] = dy
            vec_z[iter] = dz
        else
            vec_x[iter] = vec_x[iter-1] + dx
            vec_y[iter] = vec_y[iter-1] + dy
            vec_z[iter] = vec_z[iter-1] + dz
        end
    end
    return (vec_x=vec_x, vec_y=vec_y, vec_z=vec_z)
end