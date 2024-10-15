

"""
    static_mockup(parameters)

Builds a 3D mockup from a set of parameters. The mockup is static, *i.e.* it does not change over time.
"""
function static_mockup(parameters)
    rng = Random.MersenneTwister(parameters["seed"])
    # We want to avoid problem at young stage when nb_leaves_emitted < nb_fronds
    nb_fronds = length(parameters["rachis_biomass"])
    nb_internodes = max(nb_fronds + 1, parameters["nb_leaves_emitted"])
    mtg = mtg_skeleton(nb_internodes)
    mtg[1][:stem_bending] = VPalm.stem_bending(
        parameters["stem_bending_mean"],
        parameters["stem_bending_sd"]; rng=rng)
    stem_height = VPalm.stem_height(
        parameters["nb_leaves_emitted"],
        parameters["initial_stem_height"],
        parameters["stem_height_coefficient"],
        parameters["internode_length_at_maturity"],
        parameters["stem_growth_start"],
        parameters["stem_height_variation"]; rng=rng)
    stem_diameter = VPalm.stem_diameter(
        parameters["rachis_length_reference"],
        parameters["stem_diameter_max"],
        parameters["stem_diameter_slope"],
        parameters["stem_diameter_inflection"],
        parameters["stem_diameter_residual"]; rng=rng)

    for i in 1:nb_internodes
        mtg[(rank-1)*2 + 3][:Width] = VPalm.internode_diameter(
            rank,
            nb_internodes,
            stem_diameter,
            parameters["stem_base_shrinkage"],
            parameters["stem_top_shrinkage"],
            parameters["fronds_in_sheath"])
        mtg[(rank-1)*2 + 3][:Length] = VPalm.internode_length(
            rank,
            nb_internodes,
            nb_fronds,
            stem_height,
            parameters["internode_rank_no_exapansion"],
            parameters["nb_internodes_before_planting"],
            parameters["internode_final_length"])
        mtg[(rank-1)*2 + 3][:Internode_rank] = nb_internodes - i
        mtg[(rank-1)*2 + 3][:Orthotropy] = 0.05
        mtg[(rank-1)*2 + 3][:XEuler] = VPalm.phyllotactic_angle(
            parameters["phyllotactic_angle_mean"],
            parameters["phyllotactic_angle_sd"]; rng=rng)
    end
    return mtg
end