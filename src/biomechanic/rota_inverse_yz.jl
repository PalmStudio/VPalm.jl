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
    # For the inverse rotation, we use the negative angles
    # and the inverse order (ZY instead of YZ)
    # Rotations.jl convention is that RotZY means "first rotate around Y, then Z"
    # so we use RotYZ here with negated angles
    rotation = RotYZ(-agl_z, -agl_y)

    # Apply the rotation to the point
    return rotation * op
end