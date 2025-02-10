function add_geometry!(mtg, refmesh_internode, refmesh_snag)
    stem_diameter = mtg[1].stem_diameter
    stem_bending = mtg[1].stem_bending
    internode_width = stem_diameter
    snag_insertion_angle = deg2rad(-35.0) # deg2rad(-20.0 - 35.0)
    internode_height = 0.0
    snag_rotation = 0.0

    snag_width = 0.20 # see defaultOrthotropyAttribute in the trunk in the java implementation
    snag_height = 0.15
    snag_length = 3.0

    traverse!(mtg) do node
        if symbol(node) == "Internode"
            snag_rotation += deg2rad(node.XEuler)
            stem_bending += deg2rad(node.Orthotropy)
            internode_width = node.Width > 0.0 ? node.Width : 0.01
            mesh_transformation = Meshes.Scale(internode_width, internode_width, node.Length) → Meshes.Translate(0.0, 0.0, internode_height) → Meshes.Rotate(RotZ(snag_rotation)) → Meshes.Rotate(RotY(stem_bending))
            node.geometry = PlantGeom.Geometry(ref_mesh=refmesh_internode, transformation=mesh_transformation)
            internode_height += node.Length
        elseif symbol(node) == "Leaf"
            if !node.is_alive
                # Dead leaf, we keep the snag only
                mesh_transformation = Meshes.Scale(snag_length, snag_width, snag_height) → Meshes.Rotate(RotY(snag_insertion_angle)) → Meshes.Translate(internode_width, 0.0, internode_height) → Meshes.Rotate(RotZ(snag_rotation)) → Meshes.Rotate(RotY(stem_bending))
                node.geometry = PlantGeom.Geometry(ref_mesh=refmesh_snag, transformation=mesh_transformation)
            else
                nothing
            end
        elseif symbol(node) == "Petiole"
            # petiole_width = node.petiole_width
            # petiole_height = node.petiole_height
            # petiole_length = node.petiole_length
            # petiole_insertion_angle = deg2rad(node.petiole_insertion_angle)
            # petiole_rotation = deg2rad(node.petiole_rotation)

        end
    end
end