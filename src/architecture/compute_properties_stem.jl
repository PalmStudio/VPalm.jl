function compute_properties_stem!(node, parameters, rng)
    node[:stem_bending] = VPalm.stem_bending(
        parameters["stem_bending_mean"],
        parameters["stem_bending_sd"]; rng=rng
    )
    node[:stem_height] = VPalm.stem_height(
        parameters["nb_leaves_emitted"],
        parameters["initial_stem_height"],
        parameters["stem_height_coefficient"],
        parameters["internode_length_at_maturity"],
        parameters["stem_growth_start"],
        parameters["stem_height_variation"]; rng=rng
    )
    node[:stem_diameter] = VPalm.stem_diameter(
        parameters["rachis_length_reference"],
        parameters["stem_diameter_max"],
        parameters["stem_diameter_slope"],
        parameters["stem_diameter_inflection"],
        parameters["stem_diameter_residual"],
        parameters["stem_diameter_snag"]; rng=rng
    )

    return nothing
end