function rachis(petiole_node, index, scale, leaf_rank, rachis_length, height_cpoint, width_cpoint, zenithal_cpoint_angle, parameters; rng)
    rachis_node = Node(petiole_node, NodeMTG("<", "Rachis", index, scale))

    nb_segments = parameters["rachis_nb_segments"]
    points_length, points_positions, points_bending, points_deviation, points_torsion = biomechanical_properties_rachis(
        parameters["rachis_twist_initial_angle"], parameters["rachis_twist_initial_angle_sdp"],
        parameters["elastic_modulus"], parameters["shear_modulus"],
        rachis_length,
        parameters["lenflet_length_at_b_intercept"], parameters["lenflet_length_at_b_slope"], parameters["relative_position_bpoint"],
        parameters["relative_position_bpoint_sd"], parameters["relative_length_first_leaflet"], parameters["relative_length_last_leaflet"], parameters["relative_position_leaflet_max_length"],
        parameters["rachis_fresh_weigth"][leaf_rank],
        #! change the way we index in the rachis_fresh_weigth vector, because we have values for the spears too in here, so rank <= 0
        leaf_rank, height_cpoint, zenithal_cpoint_angle, nb_segments,
        parameters["height_rachis_tappering"], rng
    )

    last_parent = rachis_node
    zenithal_angle = zenithal_cpoint_angle
    azimuthal_angle = 0.0
    torsion_angle = 0.0
    for p in eachindex(points_positions)
        rachis_segment_node = Node(last_parent, NodeMTG(p == 1 ? "/" : "<", "RachisSegment", p, 6))
        rachis_segment_node.width = rachis_width(p / nb_segments, width_cpoint, parameters["rachis_width_tip"])
        rachis_segment_node.height = rachis_height(p / nb_segments, height_cpoint, parameters["height_rachis_tappering"])
        rachis_segment_node.length = points_length[p]
        rachis_segment_node.zenithal_angle = points_bending[p] - zenithal_angle
        rachis_segment_node.azimuthal_angle = points_deviation[p] - azimuthal_angle
        rachis_segment_node.torsion_angle = points_torsion[p] - torsion_angle

        zenithal_angle = rachis_segment_node.zenithal_angle
        azimuthal_angle = rachis_segment_node.azimuthal_angle
        torsion_angle = rachis_segment_node.torsion_angle
    end
end

