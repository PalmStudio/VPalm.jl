function compute_properties_leaf!(node, index, nb_internodes, nb_leaves_alive, parameters, rng)
    leaf_rank = nb_internodes - index + 1
    node[:rank] = leaf_rank

    # Is the leaf alive of dead (snag)?
    is_alive = leaf_rank <= nb_leaves_alive
    node[:is_alive] = is_alive

    if is_alive
        node[:zenithal_insertion_angle] = VPalm.leaf_insertion_angle(
            leaf_rank,
            parameters["leaf_max_angle"],
            parameters["leaf_slope_angle"],
            parameters["leaf_inflection_angle"]
        )
        #! important: we use leaf rank to index into the vector of leaf lengths but we should use another index here
        node[:rachis_length] = rachis_expansion(leaf_rank, parameters["rachis_final_lengths"][leaf_rank])
        node[:petiole_deviation_angle] = normal_deviation_draw(5.0, rng) #! make this a parameter!!!
        node[:zenithal_cpoint_angle] =
            max(
                c_point_angle(leaf_rank, parameters["cpoint_decli_intercept"], parameters["cpoint_decli_slope"], parameters["cpoint_angle_SDP"]),
                node[:zenithal_insertion_angle]
            )
        #! RV: I add this new thing were the zenithal cpoint angle cannot be lower than the insertion angle:
        # I do that because it would be weird if a leaf was going upward.
    end

    return nothing
end


"""
    rachis_expansion(leaf_rank, rachis_final_length)

    Simple function to compute the rachis expansion (using an expansion factor)
        based on the leaf rank.

    # Arguments

    - `leaf_rank`: The rank of the leaf.
    - `rachis_final_length`: The final length of the rachis.
"""
function rachis_expansion(leaf_rank, rachis_final_length)
    expansion_factor = leaf_rank < 2 ? 0.7 : 1.0
    return rachis_final_length * expansion_factor
end