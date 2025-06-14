function rachis(unique_mtg_id, index, scale, leaf_rank, rachis_length, height_cpoint, width_cpoint, zenithal_cpoint_angle, parameters; rng)
    rachis_node = Node(unique_mtg_id[], MutableNodeMTG("<", "Rachis", index, scale), Dict{Symbol,Any}())
    unique_mtg_id[] += 1

    nb_segments = parameters["rachis_nb_segments"]
    points_length, points_positions, points_bending, points_deviation, points_torsion, x, y, z = biomechanical_properties_rachis(
        parameters["rachis_twist_initial_angle"], parameters["rachis_twist_initial_angle_sdp"],
        parameters["elastic_modulus"], parameters["shear_modulus"],
        rachis_length,
        parameters["leaflet_length_at_b_intercept"], parameters["leaflet_length_at_b_slope"], parameters["relative_position_bpoint"],
        parameters["relative_position_bpoint_sd"], parameters["relative_length_first_leaflet"], parameters["relative_length_last_leaflet"], parameters["relative_position_leaflet_max_length"],
        parameters["rachis_fresh_weight"][leaf_rank], # Expected in kg
        #! change the way we index in the rachis_fresh_weight vector, because we have values for the spears too in here, so rank <= 0
        leaf_rank, height_cpoint, zenithal_cpoint_angle, nb_segments,
        parameters["height_rachis_tappering"],
        parameters["biomechanical_model"]["nb_sections"],
        parameters["biomechanical_model"]["iterations"],
        deg2rad(parameters["biomechanical_model"]["angle_max"]);
        verbose=true, rng=rng
    )

    last_parent = rachis_node
    # Computing the increment in angles between each segment (local) instead of the global angles:
    # rachis_segment_bendings = [0.0u"°"; diff(points_bending)]
    # rachis_segment_deviations = [0.0u"°"; diff(points_deviation)]
    # rachis_segment_torsions = [0.0u"°"; diff(points_torsion)]

    for p in eachindex(points_positions)
        rachis_segment_node = Node(unique_mtg_id[], last_parent, MutableNodeMTG(p == 1 ? "/" : "<", "RachisSegment", p, 6))
        unique_mtg_id[] += 1
        rachis_segment_node.width = rachis_width(p / nb_segments, width_cpoint, parameters["rachis_width_tip"])
        rachis_segment_node.height = rachis_height(p / nb_segments, height_cpoint, parameters["height_rachis_tappering"])
        rachis_segment_node.length = points_length[p]
        rachis_segment_node.zenithal_angle_global = points_bending[p]
        rachis_segment_node.azimuthal_angle_global = points_deviation[p]
        rachis_segment_node.torsion_angle_global = points_torsion[p]
        rachis_segment_node.x = x[p]
        rachis_segment_node.y = y[p]
        rachis_segment_node.z = z[p]

        last_parent = rachis_segment_node
    end

    # We force the last node to take the angles values of its parent node, because the biomechanical model can give
    # weird values at the boundaries:
    last_parent.zenithal_angle_global = parent(last_parent).zenithal_angle_global
    last_parent.azimuthal_angle_global = parent(last_parent).azimuthal_angle_global
    last_parent.torsion_angle_global = parent(last_parent).torsion_angle_global

    return rachis_node
end
