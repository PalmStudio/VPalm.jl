"""
    rota_yz(op, agl_y, agl_z)

Rotate a 3D vector around the Y and then the Z axes.

# Arguments
- `op`: A 3D vector (x, y, z).
- `agl_y`: Rotation angle around the Y axis (radians).
- `agl_z`: Rotation angle around the Z axis (radians).

# Returns
- The rotated vector.
"""
function rota_yz(op, agl_y, agl_z)
    # Rotation around OY
    cs_y = cos(agl_y)
    sn_y = sin(agl_y)

    mat_rot_y = [cs_y 0 -sn_y;
        0 1 0;
        sn_y 0 cs_y]

    vec_rot_y = mat_rot_y * op

    # Rotation around OZ
    cs_z = cos(agl_z)
    sn_z = sin(agl_z)

    mat_rot_z = [cs_z -sn_z 0;
        sn_z cs_z 0;
        0 0 1]

    return mat_rot_z * vec_rot_y
end
#! We should use the Rotations package instead