function add_geometry!(mtg, refmesh_cylinder, refmesh_snag)
    stem_diameter = mtg[1].stem_diameter
    stem_bending = mtg[1].stem_bending
    internode_width = stem_diameter
    snag_insertion_angle = deg2rad(-35.0) # deg2rad(-20.0 - 35.0)
    internode_height = 0.0u"m"
    snag_rotation = 0.0

    snag_width = 0.20u"m" # see defaultOrthotropyAttribute in the trunk in the java implementation
    snag_height = 0.15u"m"
    snag_length = 3.0u"m"
    position_section = Ref(Meshes.Point(0.0, 0.0, 0.0))
    angles = [0.0, 0.0, 0.0]

    traverse!(mtg, symbol=["Internode", "Leaf", "Petiole", "Rachis"]) do node
        if symbol(node) == "Internode"
            snag_rotation += deg2rad(node.XEuler)
            stem_bending += deg2rad(node.Orthotropy)
            internode_width = node.Width > 0.0u"m" ? node.Width : 0.01u"m"
            mesh_transformation = Meshes.Scale(ustrip(internode_width), ustrip(internode_width), ustrip(node.Length)) → Meshes.Translate(0.0u"m", 0.0u"m", internode_height) → Meshes.Rotate(RotZ(snag_rotation)) → Meshes.Rotate(RotY(stem_bending)) #! rotations should be merged
            node.geometry = PlantGeom.Geometry(ref_mesh=refmesh_cylinder, transformation=mesh_transformation)
            internode_height += node.Length
        elseif symbol(node) == "Leaf"
            if !node.is_alive
                # Dead leaf, we keep the snag only
                mesh_transformation = Meshes.Scale(ustrip(snag_length), ustrip(snag_width), ustrip(snag_height)) → Meshes.Rotate(RotY(snag_insertion_angle)) → Meshes.Translate(internode_width, 0.0u"m", internode_height) → Meshes.Rotate(RotZ(snag_rotation)) → Meshes.Rotate(RotY(stem_bending))
                node.geometry = PlantGeom.Geometry(ref_mesh=refmesh_snag, transformation=mesh_transformation)
            else
                nothing
            end
        elseif symbol(node) == "Petiole"
            # Initialise the position and angles for the petiole to 0.0
            position_section[] = Meshes.Point(0.0, 0.0, 0.0)
            angles .= 0.0
            add_section_geometry!(node, internode_width, internode_height, snag_rotation, stem_bending, refmesh_cylinder, "PetioleSegment", position_section, angles)
        elseif symbol(node) == "Rachis"
            add_section_geometry!(node, internode_width, internode_height, snag_rotation, stem_bending, refmesh_cylinder, "RachisSegment", position_section, angles)
            # Note: we use the position and angles of the last petiole section to initialize the rachis
        end
    end
end