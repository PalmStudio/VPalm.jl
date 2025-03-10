"""
    mtg_skeleton(nb_internodes)

Makes an MTG skeleton with `nb_leaves_emitted` leaves, including all intermediate organs:

- Plant: the whole palm
- Stem: the stem of the plant, *i.e.* the remaining part of the plant after the leaves have been removed
- Phytomer: the part that includes the leaf and the internode
- Internodes: the part of the phytomer that is between two leaves
- Leaf: the leaf of the plant, also called frond

Note: this skeleton does not include reproductive organs (inflorescences, fruits) or the scales that decompose the leaf (petiole, rachis, leaflets).

# Arguments

- `nb_internodes`: The number of internodes to emit.

# Examples

```julia
file = joinpath(dirname(dirname(pathof(VPalm))), "test", "files", "parameter_file.yml")
parameters = read_parameters(file)
mtg_skeleton(parameters)
```
"""
function mtg_skeleton(parameters; rng=Random.MersenneTwister(parameters["seed"]))
    nb_internodes = parameters["nb_leaves_emitted"] + parameters["nb_internodes_before_planting"] # The number of internodes emitted since the seed
    nb_leaves_alive = floor(Int, mean_and_sd(parameters["nb_leaves_mean"], parameters["nb_leaves_sd"]; rng=rng))
    nb_leaves_alive = min(nb_leaves_alive, nb_internodes)
    nb_petiole_segments = parameters["petiole_nb_segments"]

    @assert length(parameters["rachis_fresh_weight"]) >= nb_leaves_alive "The number of rachis biomass values should be greater than or equal to the number of leaves alive ($nb_leaves_alive)."

    # Plant / Scale 1
    plant = Node(NodeMTG("/", "Plant", 1, 1))

    # Stem (& Roots) / Scale 2
    #roots = Node(plant, NodeMTG("+", "RootSystem", 1, 2))
    stem = Node(plant, NodeMTG("+", "Stem", 1, 2))
    compute_properties_stem!(stem, parameters, rng)

    stem_height = stem[:stem_height]
    stem_diameter = stem[:stem_diameter]

    # Phytomer / Scale 3
    phytomer = Node(stem, NodeMTG("/", "Phytomer", 1, 3))

    # Internode & Leaf / Scale 4
    internode = Node(phytomer, NodeMTG("/", "Internode", 1, 4))
    compute_properties_internode!(internode, 1, nb_internodes, nb_leaves_alive, stem_height, stem_diameter, parameters, rng)
    leaf = Node(internode, NodeMTG("+", "Leaf", 1, 4))
    compute_properties_leaf!(leaf, 1, nb_internodes, nb_leaves_alive, parameters, rng)

    # Loop on internodes
    for i in 2:nb_internodes
        phytomer = Node(phytomer, NodeMTG("<", "Phytomer", i, 3))
        internode = Node(phytomer, NodeMTG("/", "Internode", i, 4))
        compute_properties_internode!(internode, i, nb_internodes, nb_leaves_alive, stem_height, stem_diameter, parameters, rng)
        leaf = Node(internode, NodeMTG("+", "Leaf", i, 4))
        compute_properties_leaf!(leaf, i, nb_internodes, nb_leaves_alive, parameters, rng)
        # Loop on present leaves
        if leaf[:is_alive]
            # Build the petiole
            petiole_node = petiole(leaf, i, 5, leaf.rachis_length, leaf.zenithal_insertion_angle, leaf.zenithal_cpoint_angle, parameters; rng=rng)
            # Build the rachis
            rachis_node = rachis(petiole_node, i, 5, leaf.rank, leaf.rachis_length, petiole_node.height_cpoint, petiole_node.width_cpoint, leaf.zenithal_cpoint_angle, parameters; rng=rng)
            # Add the leaflets to the rachis:
            leaflets(rachis_node, i, 5, leaf.rank, leaf.rachis_length, petiole_node.height_cpoint, petiole_node.width_cpoint, leaf.zenithal_cpoint_angle, parameters; rng=rng)
        end

        # add leaflets, ls
    end

    # Compute the geometry of the plant
    # Note: we could do this at the same time than the architecture, but it is separated here for clarity. The downside is that we traverse the mtg twice, but it is pretty cheap.
    refmesh_cylinder = PlantGeom.RefMesh("cylinder", VPalm.cylinder())
    refmesh_snag = PlantGeom.RefMesh("Snag", VPalm.snag(0.05, 1.0, 1.0))

    add_geometry!(plant, refmesh_cylinder, refmesh_snag)

    return plant
end
