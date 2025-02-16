function rachis(petiole_node, index, scale, leaf_rank, rachis_length, height_cpoint, width_cpoint, zenithal_cpoint_angle, parameters; rng)
    rachis_node = Node(petiole_node, NodeMTG("<", "Rachis", index, scale))

    nb_segments = parameters["rachis_nb_segments"]
    points_length, points_positions, points_bending, points_deviation, points_torsion, x, y, z = biomechanical_properties_rachis(
        parameters["rachis_twist_initial_angle"], parameters["rachis_twist_initial_angle_sdp"],
        parameters["elastic_modulus"], parameters["shear_modulus"],
        rachis_length,
        parameters["lenflet_length_at_b_intercept"], parameters["leaflet_length_at_b_slope"], parameters["relative_position_bpoint"],
        parameters["relative_position_bpoint_sd"], parameters["relative_length_first_leaflet"], parameters["relative_length_last_leaflet"], parameters["relative_position_leaflet_max_length"],
        parameters["rachis_fresh_weigth"][leaf_rank] / 1000.0, # Expected in kg
        #! change the way we index in the rachis_fresh_weigth vector, because we have values for the spears too in here, so rank <= 0
        leaf_rank, height_cpoint, zenithal_cpoint_angle, nb_segments,
        parameters["height_rachis_tappering"], rng
    )

    last_parent = rachis_node
    zenithal_angle_prev = -zenithal_cpoint_angle
    azimuthal_angle_prev = 0.0
    torsion_angle_prev = 0.0
    for p in eachindex(points_positions)
        rachis_segment_node = Node(last_parent, NodeMTG(p == 1 ? "/" : "<", "RachisSegment", p, 6))
        rachis_segment_node.width = rachis_width(p / nb_segments, width_cpoint, parameters["rachis_width_tip"])
        rachis_segment_node.height = rachis_height(p / nb_segments, height_cpoint, parameters["height_rachis_tappering"])
        rachis_segment_node.length = points_length[p]

        rachis_segment_node.zenithal_angle = points_bending[p] - zenithal_angle_prev
        rachis_segment_node.zenithal_angle_bending = points_bending[p]
        zenithal_angle_prev = points_bending[p]
        rachis_segment_node.azimuthal_angle = points_deviation[p] - azimuthal_angle_prev
        azimuthal_angle_prev = points_deviation[p]
        rachis_segment_node.torsion_angle = points_torsion[p] - torsion_angle_prev
        torsion_angle_prev = points_torsion[p]
        rachis_segment_node.x = x[p]
        rachis_segment_node.y = y[p]
        rachis_segment_node.z = z[p]

        last_parent = rachis_segment_node
    end

    # We force the last node to take the angles values of its parent node, because the biomechanical model can give
    # weird values at the boundaries:
    last_parent.zenithal_angle = parent(last_parent).zenithal_angle
    last_parent.azimuthal_angle = parent(last_parent).azimuthal_angle
    last_parent.torsion_angle = parent(last_parent).torsion_angle

end

