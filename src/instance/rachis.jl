function rachis(petiole_node, index, scale, leaf_rank, rachis_length, zenithal_cpoint_angle, parameters; rng=rng)
    rachis_node = Node(petiole_node, NodeMTG("<", "Rachis", index, scale))

    compute_properties_rachis!(
        rachis_node,
        index,
        parameters["rachis_twist_initial_angle"], parameters["rachis_twist_initial_angle_sdp"],
        parameters["elastic_modulus"], parameters["shear_modulus"],
        rachis_length,
        parameters["lenflet_length_at_b_intercept"], parameters["lenflet_length_at_b_slope"], parameters["relative_position_bpoint"],
        parameters["relative_position_bpoint_sd"], parameters["relative_length_first_leaflet"], parameters["relative_length_last_leaflet"], parameters["relative_position_leaflet_max_length"],
        parameters["rachis_fresh_weigth"], leaf_rank, petiole_node.height_cpoint, zenithal_cpoint_angle, parameters["rachis_nb_segments"],
        parameters["height_rachis_tappering"];
        rng=rng
    )

end