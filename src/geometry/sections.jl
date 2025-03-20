"""
    add_section_geometry!(
        node, internode_width, internode_height, internode_phyllotaxy, stem_bending, 
        refmesh_cylinder, position_section=Ref(Meshes.Point(0.0, 0.0, 0.0)), angles=[0.0, 0.0, 0.0],
        type::String,
    )

Create the petiole/rachis sections geometry based on their dimensions.

# Arguments

- `node`: the MTG node of the petiole/rachis
- `refmesh_cylinder`: the reference mesh used for a cylinder (`PlantGeom.RefMesh`)
- `internode_width`: the width of the internode on the stipe (m)
- `internode_height`: the height of the internode on the stipe (m)
- `internode_phyllotaxy`: the phyllotaxy of the internode on the stipe (°)
- `stem_bending`: the bending of the stipe (°)
- `type::String`: the type of the section (`"PetioleSegment"` or `"RachisSegment"`)
- `position_section=Ref(Meshes.Point(0.0, 0.0, 0.0))`: the position of the section relative to the first one.
- `angles=[0.0u"°", 0.0u"°", 0.0u"°"]`: the angles of the section relative to the first one.
"""
function add_section_geometry!(
    node, refmesh_cylinder, internode_width=0.0u"m", internode_height=0.0u"m", internode_phyllotaxy=0.0u"°", stem_bending=0.0u"°",
    type=nothing, position_section=Ref(Meshes.Point(0.0, 0.0, 0.0)), angles=[0.0u"°", 0.0u"°", 0.0u"°"],
)

    # Check units and convert to meters and degrees:
    internode_width = @check_unit internode_width u"m"
    internode_height = @check_unit internode_height u"m"
    internode_phyllotaxy = @check_unit internode_phyllotaxy u"°"
    stem_bending = @check_unit stem_bending u"°"
    angles[1] = @check_unit angles[1] u"°"
    angles[2] = @check_unit angles[2] u"°"
    angles[3] = @check_unit angles[3] u"°"

    # If type is not provided, we use the symbol of the first child node to filter:
    type = type === nothing ? symbol(node[1]) : type
    traverse!(node, symbol=type) do node_section
        angles[1] += node_section.zenithal_angle
        angles[2] += node_section.azimuthal_angle
        angles[3] += node_section.torsion_angle
        α = angles[1] # Elevation (i.e. zenithal angle, in degrees)
        β = angles[2] # Azimuth (in degrees)
        γ = angles[3] # Torsion (in degrees)

        # Rotation matrix for the section
        rotation = RotXYZ(-deg2rad(α), deg2rad(β), deg2rad(γ))
        # rotation = RotXYZ(-deg2rad(α), deg2rad(γ), deg2rad(β))

        mesh_transformation =
            Meshes.Scale(ustrip(node_section.width) / 2.0, ustrip(node_section.height) / 2.0, ustrip(node_section.length)) →
            Meshes.Rotate(rotation) →
            Meshes.Translate(Meshes.to(position_section[])...) →
            Meshes.Rotate(RotZ(-π / 2)) → # orient the reference cylinder to face X forward
            # Positioning along the stem:
            Meshes.Translate(internode_width, zero(internode_width), internode_height) →
            Meshes.Rotate(RotZ(internode_phyllotaxy)) →
            Meshes.Rotate(RotY(stem_bending))

        node_section.geometry = PlantGeom.Geometry(ref_mesh=refmesh_cylinder, transformation=mesh_transformation)

        x = cosd(α) * cosd(β)
        y = sind(α) * cosd(β)
        z = sind(β)
        n = Meshes.Vec(z, y, x) # vector normal to the plane at point p1

        position_section[] += n * ustrip(node_section.length)
    end

    return nothing
end
