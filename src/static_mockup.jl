

"""
    static_mockup(parameters)

Builds a 3D mockup from a set of parameters. The mockup is static, *i.e.* it does not change over time.
"""
function static_mockup(parameters)
    mockup = mtg_skeleton(parameters["nb_leaves_emitted"])
end