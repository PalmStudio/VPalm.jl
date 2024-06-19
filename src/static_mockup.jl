

"""
    static_mockup(parameters)

Builds a 3D mockup from a set of parameters. The mockup is static, *i.e.* it does not change over time.
"""
function static_mockup(parameters)
    rng = Random.MersenneTwister(parameters["seed"])
    # To avoid problem at young stage when nb_leaves_emitted < nb_fronds, we define nb_internodes
    # This comes from Vpalm java implementation
    # Note: why not just use max(nb_fronds, nb_leaves_emitted) instead of max(nb_fronds + 1, nb_leaves_emitted)?
    nb_internodes = max(parameters["nb_fronds"] +1 , parameters["nb_leaves_emitted"])
    mtg = mtg_skeleton(parameters["nb_internodes"])
    mtg[1][:stem_bending] = stem_bending(
        parameters["stem_bending_mean"],
        parameters["stem_bending_sd"]
        )
    mtg[1][:stem_height] = stem_height(
        parameters["nb_leaves_emitted"],
        parameters["initial_stem_height"],
        parameters["stem_height_coefficient"],
        parameters["internode_length_at_maturity"],
        parameters["stem_growth_start"],
        parameters["stem_height_variation"]
        )
    mtg[1][:stem_diameter] = stem_diameter(
        parameters["rachis_reference_length"],
        parameters["stem_diameter_max"],
        parameters["stem_diameter_slope"],
        parameters["stem_diameter_inflection"],
        parameters["stem_diameter_variation"]
        )
    for i in 1: nb_leaves_emitted

    end
end