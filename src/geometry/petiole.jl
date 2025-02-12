"""
    add_petiole_section_geometry!(petiole_node, internode_width, internode_height, refmesh_cylinder)

Create the petiole sections geometry based on their dimensions.

# Arguments

- `petiole_node`: the MTG node of the petiole
- `internode_width`: the width of the internode on the stipe
- `internode_height`: the heigth of the internode on the stipe
- `refmesh_cylinder`: the reference mesh used for a cylinder
"""
function add_petiole_section_geometry!(petiole_node, internode_width, internode_height, refmesh_cylinder)
    section_insertion_angle = 0.0
    section_azimuthal_angle = 0.0
    position_section = Ref(Meshes.Point(internode_width, internode_height, 0.0))
    # top_point = Point(0.0, 0.0, 0.0)
    traverse!(petiole_node[1]) do node_section
        section_insertion_angle += node_section.zenithal_angle
        section_azimuthal_angle += node_section.azimuthal_angle

        # section_dimensions = [node_section.length, node_section.width, node_section.height]
        # section_dimensions = [node_section.width, node_section.length, node_section.height]
        # section_dimensions = [node_section.length, node_section.height, node_section.width]

        # section_dimensions = [node_section.height, node_section.width, node_section.length]
        section_dimensions = [node_section.width, node_section.height, node_section.length]

        mesh_transformation =
            Meshes.Scale(section_dimensions...) →
            Meshes.Rotate(RotX(π / 2 + deg2rad(section_insertion_angle))) → # We first rotate the reference mesh to get X = width, Y = length, Z = height
            Meshes.Rotate(Rotations.RotY(deg2rad(section_azimuthal_angle))) →
            # Meshes.Rotate(Rotations.RotX(deg2rad(section_insertion_angle))) →
            Meshes.Translate(Meshes.to(position_section[])...)
        # Meshes.Translate(0.0, 0.0, 1.3)

        # Meshes.Rotate(RotZ(snag_rotation)) → Meshes.Rotate(RotY(stem_bending))
        node_section.geometry = PlantGeom.Geometry(ref_mesh=refmesh_cylinder, transformation=mesh_transformation)

        # Moving the next refmesh by the lenght of the previous mesh
        # position_section[2] += node_section.length

        α = section_insertion_angle # Elevation (i.e. zenithal angle, in degrees)
        β = section_azimuthal_angle # Azimuth (in degrees)
        x = cosd(α) * cosd(β)
        y = sind(α) * cosd(β)
        z = sind(β)
        # n = Meshes.Vec(x, y, z) # vector normal to the plane at point p1
        # n = Meshes.Vec(x, z, y) # vector normal to the plane at point p1
        # n = Meshes.Vec(z, y, x) # vector normal to the plane at point p1
        n = Meshes.Vec(z, x, y) # vector normal to the plane at point p1
        position_section[] += n * node_section.length
        # position_section[] = position_section[] + n * node_section.length #- Meshes.Vec(0.0, 0.0, node_section.width / 2)
        # position_section[] = position_section[] + Meshes.Vec(0.0, node_section.length, node_section.width / 2)
    end
end
