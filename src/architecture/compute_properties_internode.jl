function compute_properties_internode!(node, index, nb_internodes, nb_leaves_alive, stem_height, stem_diameter, parameters, rng)
    node[:Width] = VPalm.internode_diameter(
        index,
        nb_internodes,
        stem_diameter,
        parameters["stem_base_shrinkage"],
        parameters["stem_top_shrinkage"],
        parameters["nb_leaves_in_sheath"]
    )

    node[:Length] = internode_length(
        index,
        nb_internodes,
        stem_height,
        parameters["internode_rank_no_expansion"],
        parameters["nb_internodes_before_planting"],
        0.001
        # internode_min_height
    )

    node[:rank] = nb_internodes - index + 1
    node[:Orthotropy] = 0.05
    node[:XEuler] = VPalm.phyllotactic_angle(
        parameters["phyllotactic_angle_mean"],
        parameters["phyllotactic_angle_sd"]; rng=rng
    )
    return nothing
end