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
    # Use the combined RotYZ rotation directly
    rotation = RotYZ(agl_y, agl_z)

    # Apply the rotation to the point
    return rotation * op
end