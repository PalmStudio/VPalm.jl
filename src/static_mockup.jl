

"""
    static_mockup(parameters)

Builds a 3D mockup from a set of parameters. The mockup is static, *i.e.* it does not change over time.
"""
function static_mockup(parameters)
    rng = Random.MersenneTwister(parameters["seed"])
    nb_leaves = length(parameters["rachis_biomass"])
    # We want to avoid problem at young stage when nb_leaves_emitted < nb_leaves
    nb_internodes = max(nb_leaves + 1, parameters["nb_leaves_emitted"])
    mtg = mtg_skeleton(nb_internodes)
    mtg[1][:stem_bending] = VPalm.stem_bending(
        parameters["stem_bending_mean"],
        parameters["stem_bending_sd"]; rng=rng
    )
    stem_height = VPalm.stem_height(
        parameters["nb_leaves_emitted"],
        parameters["initial_stem_height"],
        parameters["stem_height_coefficient"],
        parameters["internode_length_at_maturity"],
        parameters["stem_growth_start"],
        parameters["stem_height_variation"]; rng=rng
    )
    stem_diameter = VPalm.stem_diameter(
        parameters["rachis_length_reference"],
        parameters["stem_diameter_max"],
        parameters["stem_diameter_slope"],
        parameters["stem_diameter_inflection"],
        parameters["stem_diameter_residual"],
        parameters["stem_diameter_snag"]; rng=rng
    )

    for rank in 1:nb_internodes
        # First retrieve the Internode and Leaf Nodes
        internode = get_node(mtg,(rank-1)*3 + 4)
        leaf = get_node(mtg,(rank-1)*3 + 5)
        # Then associate attributes
        internode[:Width] = VPalm.internode_diameter(
            rank,
            nb_internodes,
            stem_diameter,
            parameters["stem_base_shrinkage"],
            parameters["stem_top_shrinkage"],
            parameters["nb_leaves_in_sheath"]
        )
        internode[:Length] = VPalm.internode_length(
            rank,
            nb_internodes,
            nb_leaves,
            stem_height,
            parameters["internode_rank_no_expansion"],
            parameters["nb_internodes_before_planting"],
            parameters["internode_final_length"]
        )
        leaf_rank = nb_internodes - rank
        internode[:rank] = leaf_rank
        internode[:Orthotropy] = 0.05
        internode[:XEuler] = VPalm.phyllotactic_angle(
            parameters["phyllotactic_angle_mean"],
            parameters["phyllotactic_angle_sd"]; rng=rng
        )
        leaf[:rank] = leaf_rank
        leaf[:rachis_freshweight] = parameters["rachis_biomass"][rank]
        leaf[:YInsertionAngle] = VPalm.leaf_insertion_angle(
            leaf_rank,
            parameters["leaf_max_angle"],
            parameters["leaf_slope_angle"],
            parameters["leaf_inflection_angle"]
        )
    end
    return mtg
end