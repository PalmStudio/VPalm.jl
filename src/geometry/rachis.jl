"""
    add_section_geometry!(
        node, internode_width, internode_height, internode_phyllotaxy, stem_bending, 
        refmesh_cylinder, position_section=Ref(Meshes.Point(0.0, 0.0, 0.0)
    )

Create the petiole/rachis sections geometry based on their dimensions.

# Arguments

- `node`: the MTG node of the petiole/rachis
- `internode_width`: the width of the internode on the stipe (m)
- `internode_height`: the heigth of the internode on the stipe (m)
- `internode_phyllotaxy`: the phyllotaxy of the internode on the stipe (°)
- `stem_bending`: the bending of the stipe (°)
- `refmesh_cylinder`: the reference mesh used for a cylinder (`PlantGeom.RefMesh`)
- `position_section=Ref(Meshes.Point(0.0, 0.0, 0.0)`: the position of the section relative to the first one.
"""
function add_section_geometry!(node, internode_width, internode_height, internode_phyllotaxy, stem_bending, refmesh_cylinder, position_section=Ref(Meshes.Point(0.0, 0.0, 0.0))
)
    section_insertion_angle = 0.0
    section_azimuthal_angle = 0.0
    section_twist_angle = 0.0

    traverse!(node[1]) do node_section
        section_insertion_angle += node_section.zenithal_angle
        section_azimuthal_angle += node_section.azimuthal_angle
        section_twist_angle += node_section.torsion_angle

        section_dimensions = [node_section.width / 2.0, node_section.height / 2.0, node_section.length]

        mesh_transformation =
            Meshes.Scale(section_dimensions...) →
            Meshes.Rotate(RotXYZ(-deg2rad(section_insertion_angle), deg2rad(section_azimuthal_angle), deg2rad(section_twist_angle))) →
            Meshes.Translate(Meshes.to(position_section[])...) →
            Meshes.Rotate(RotZ(-π / 2)) → # orient the reference cylinder to face X forward
            # Positioning along the stem:
            Meshes.Translate(internode_width, 0.0, internode_height) →
            Meshes.Rotate(RotZ(internode_phyllotaxy)) →
            Meshes.Rotate(RotY(stem_bending))

        # Meshes.Rotate(RotZ(snag_rotation)) → Meshes.Rotate(RotY(stem_bending))
        node_section.geometry = PlantGeom.Geometry(ref_mesh=refmesh_cylinder, transformation=mesh_transformation)

        α = section_insertion_angle # Elevation (i.e. zenithal angle, in degrees)
        β = section_azimuthal_angle # Azimuth (in degrees)
        x = cosd(α) * cosd(β)
        y = sind(α) * cosd(β)
        z = sind(β)
        n = Meshes.Vec(z, y, x) # vector normal to the plane at point p1
        position_section[] += n * node_section.length
    end

    return nothing
end
