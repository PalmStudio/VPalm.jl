"""
    agl_vers_xyz(dist_p2p1, vangle_xy, vangle_xz)

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
function agl_vers_xyz(dist_p2p1, vangle_xy, vangle_xz)
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