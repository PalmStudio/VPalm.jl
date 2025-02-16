"""
    xyz_vers_agl(vec_x, vec_y, vec_z)

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
function xyz_vers_agl(vec_x, vec_y, vec_z)
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

        if p2p1[1] != 0
            # Distances
            dist_p2p1[iter] = sqrt(p2p1[1]^2 + p2p1[2]^2 + p2p1[3]^2)

            # Angles
            vangle_xy[iter] = atan(p2p1[3] / sqrt(p2p1[1]^2 + p2p1[2]^2))
            vangle_xz[iter] = atan(p2p1[2] / p2p1[1])
        end
    end

    return (dist_p2p1=dist_p2p1, vangle_xy=vangle_xy, vangle_xz=vangle_xz)
end