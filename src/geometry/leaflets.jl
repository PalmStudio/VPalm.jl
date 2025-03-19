"""
    add_leaflet_geometry!(
        leaflet_node, 
        rachis_position, 
        rachis_orientation, 
        rachis_rotation, 
        stem_bending, 
        refmesh_plane
    )

Create the leaflet geometry based on its segments.

# Arguments
- `leaflet_node`: The MTG node of the leaflet
- `rachis_position`: Position of the rachis section where the leaflet is attached
- `rachis_orientation`: Orientation angles [zenithal, azimuthal, torsion] of the rachis section
- `rachis_rotation`: Rotation of the rachis due to phyllotaxy (radians)
- `stem_bending`: Bending of the stem (radians)
- `refmesh_plane`: Reference mesh used for the planar leaflet segments

# Returns
- Nothing (the geometry is added directly to the leaflet node and its segments)
"""
function add_leaflet_geometry!(
    leaflet_node,
    rachis_position,
    rachis_orientation,
    rachis_rotation,
    stem_bending,
    refmesh_plane
)
    # Extract basic leaflet properties
    side = leaflet_node["side"]
    h_angle = deg2rad(leaflet_node["azimuthal_angle"]) # Horizontal angle (insertion angle in Z)  
    v_angle = deg2rad(leaflet_node["zenithal_angle"]) # Vertical angle (insertion angle in X)
    torsion = deg2rad(leaflet_node["torsion_angle"])    # Twist around leaflet's axis

    # Create reference point for the leaflet base
    leaflet_base = rachis_position

    # Accumulate position and angle for segments
    position = Ref(Meshes.Point(0.0, 0.0, 0.0))

    # Calculate the orientation for the full leaflet
    # 1. Apply leaflet's insertion angles (horizontal and vertical)
    # 2. Apply torsion around leaflet's own axis
    # 3. Apply rachis orientation (inherit from parent)
    leaflet_orientation = [
        rachis_orientation[1] + v_angle, # Add vertical insertion angle to rachis zenithal angle
        rachis_orientation[2],           # Keep rachis azimuthal angle
        rachis_orientation[3] + h_angle  # Add horizontal insertion angle to rachis rotation
    ]

    # Create a reference for the previous segment angle for segment-to-segment bending
    previous_segment_angle = 0.0

    # Process each leaflet segment
    traverse!(leaflet_node, symbol="LeafletSegment") do segment
        # Get segment properties
        width = segment["width"]
        length = segment["length"]

        # Calculate the lamina angle (the V-shape of the leaflet)
        lamina_angle = deg2rad(30.0)  # Default 30-degree angle as in Java (from ElaeisArchiTree.java)

        # Apply stiffness angle if available (segment bending due to weight)
        segment_stiffness_angle = 0.0
        if haskey(segment, "stiffness_angle")
            segment_stiffness_angle = deg2rad(segment["stiffness_angle"])
        end

        # Calculate the absolute angle of this segment by adding the stiffness angle to previous segment angle
        segment_angle = previous_segment_angle + segment_stiffness_angle
        previous_segment_angle = segment_angle

        # Calculate transformation for this segment
        # Based on ElaeisArchiTree.java line 246-254, the leaflet shape is a V-shaped plane
        mesh_transformation =
        # Scale to correct dimensions
            Meshes.Scale(ustrip(width), ustrip(width) * sin(lamina_angle / 2), ustrip(length)) →
            # Apply segment bending angle
            Meshes.Rotate(RotX(segment_angle)) →
            # Apply position offset from previous segments
            Meshes.Translate(Meshes.to(position[])...) →
            # Apply leaflet insertion and orientation
            Meshes.Rotate(RotXYZ(leaflet_orientation[1], leaflet_orientation[2], leaflet_orientation[3])) →
            # Apply torsion around leaflet axis
            Meshes.Rotate(RotX(torsion * side)) →
            # Translate to rachis position
            Meshes.Translate(Meshes.to(leaflet_base)...) →
            # Apply stem global transforms
            Meshes.Rotate(RotZ(rachis_rotation)) →  #! do we really need this? This is the phyllotaxy, but we should already have it in the rachis segment orientation
            Meshes.Rotate(RotY(stem_bending)) #! do we really need this?

        # Assign geometry to the segment
        segment.geometry = PlantGeom.Geometry(ref_mesh=refmesh_plane, transformation=mesh_transformation)

        # Update position for next segment by moving along the segment's length
        # Considering the segment's bending angle
        dx = length * cos(segment_angle)
        dy = 0.0  # No sideways movement along segment's own axis
        dz = length * sin(segment_angle)

        # Update position reference for the next segment
        position[] += Meshes.Vec(dx, dy, dz)
    end

    return nothing
end
