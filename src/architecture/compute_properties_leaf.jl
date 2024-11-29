function compute_properties_leaf!(node, index, nb_internodes, nb_leaves_alive, parameters, rng)
    leaf_rank = nb_internodes - index + 1
    node[:rank] = leaf_rank

    # Is the leaf alive of dead (snag)?
    is_alive = leaf_rank <= nb_leaves_alive
    node[:is_alive] = is_alive

    if is_alive
        node[:rachis_freshweight] = parameters["rachis_biomass"][leaf_rank]
        node[:YInsertionAngle] = VPalm.leaf_insertion_angle(
            leaf_rank,
            parameters["leaf_max_angle"],
            parameters["leaf_slope_angle"],
            parameters["leaf_inflection_angle"]
        )
    end

    return nothing
end