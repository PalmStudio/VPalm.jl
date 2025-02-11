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
    position_section = [internode_width, 0.0, internode_height]
    traverse!(petiole_node[1]) do node_section
        section_insertion_angle += node_section.zenithal_angle
        section_azimuthal_angle += node_section.azimuthal_angle

        # section_dimensions = [node_section.length, node_section.width, node_section.height]
        # section_dimensions = [node_section.width, node_section.length, node_section.height]
        # section_dimensions = [node_section.length, node_section.height, node_section.width]

        # section_dimensions = [node_section.height, node_section.width, node_section.length]
        section_dimensions = [node_section.width, node_section.length, node_section.height]

        mesh_transformation =
            Meshes.Rotate(RotX(π / 2)) → # We first rotate the reference mesh to get X = width, Y = length, Z = height
            Meshes.Scale(section_dimensions...) →
            Meshes.Translate(position_section...) →
            Meshes.Rotate(RotXY(deg2rad(section_insertion_angle), deg2rad(section_azimuthal_angle)))

        # Meshes.Rotate(RotZ(snag_rotation)) → Meshes.Rotate(RotY(stem_bending))
        node_section.geometry = PlantGeom.Geometry(ref_mesh=refmesh_cylinder, transformation=mesh_transformation)

        # Moving the next refmesh by the lenght of the previous mesh
        position_section[2] += node_section.length
    end
end
