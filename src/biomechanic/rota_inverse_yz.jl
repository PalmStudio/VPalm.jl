"""
    rota_inverse_yz(op, agl_y, agl_z)

Rotate a point around the Z axis and then the Y axis (inverse rotation).

# Arguments
- `op`: A 3D point (x, y, z) as a Vector.
- `agl_y`: Rotation angle around the Y axis (radians).
- `agl_z`: Rotation angle around the Z axis (radians).

# Returns
- The rotated point as a Vector.
"""
function rota_inverse_yz(op, agl_y, agl_z)
    agl_y = -agl_y
    agl_z = -agl_z

    # Rotation around OZ
    cs_z = cos(agl_z)
    sn_z = sin(agl_z)

    mat_rot_z = [cs_z -sn_z 0;
        sn_z cs_z 0;
        0 0 1]

    vec_rot_z = mat_rot_z * op

    # Rotation around OY
    cs_y = cos(agl_y)
    sn_y = sin(agl_y)

    mat_rot_y = [cs_y 0 -sn_y;
        0 1 0;
        sn_y 0 cs_y]

    return mat_rot_y * vec_rot_z
end