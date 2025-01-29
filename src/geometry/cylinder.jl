"""
    cylinder()
    cylinder(r, l)

Returns a normalized cylinder mesh, or a cylinder with radius `r` and length `l`.

# Arguments

- `r`: The radius of the cylinder.
- `l`: The length of the cylinder.
"""
cylinder() = Meshes.CylinderSurface(1.0) |> Meshes.discretize |> Meshes.simplexify
cylinder(r, l) = Meshes.CylinderSurface(Meshes.Point(0.0, 0.0, 0.0), Meshes.Point(0.0, 0.0, l), r) |> Meshes.discretize |> Meshes.simplexify