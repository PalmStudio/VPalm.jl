

"""
    static_mockup(parameters)

Builds a 3D mockup from a set of parameters. The mockup is static, *i.e.* it does not change over time.
"""
function static_mockup(parameters)
    rng = Random.MersenneTwister(parameters["seed"])
    mtg = mtg_skeleton(parameters["nb_leaves_emitted"])
    mtg[1][:stem_bending] = VPalm.stem_bending(
        parameters["stem_bending_mean"],
        parameters["stem_bending_sd"]; rng=rng)
    mtg[1][:stem_height] = VPalm.stem_height(
        parameters["nb_leaves_emitted"],
        parameters["initial_stem_height"],
        parameters["stem_height_coefficient"],
        parameters["internode_length_at_maturity"],
        parameters["stem_growth_start"],
        parameters["stem_height_variation"]; rng=rng)

    return mtg
end