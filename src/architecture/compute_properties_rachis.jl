function compute_properties_rachis!(
    leaf_node,
    rachis_node,
    index,
    rachis_final_length,
    c_point_width_intercept,
    c_point_width_slope,
    rng)

    leaf_node[:rachis_length] = rachis_expansion(leaf_rank, rachis_final_length)
    rachis_node[:deviation_angle] = normalDeviationDraw(5, rng)

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

