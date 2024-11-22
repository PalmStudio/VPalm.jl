"""
    mtg_skeleton(nb_leaves_emitted)

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
mtg_skeleton(3)
```
"""
function mtg_skeleton(nb_internodes, parameters; rng=Random.MersenneTwister(parameters["seed"]))
    nb_internodes = parameters["nb_leaves_emitted"]
    nb_leaves_alive = mean_and_sd(parameters["nb_leaves_mean"], parameters["nb_leaves_sd"]; rng=rng)
    nb_leaves_alive = min(nb_leaves_alive, nb_internodes)
    nb_petiole_segments = parameters["petiole_nb_segments"]
    nb_rachis_segments = parameters["rachis_nb_segments"]

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
            # Petiole / Scale 5
            petiole = Node(leaf, NodeMTG("/", "Petiole", i, 5))
            compute_properties_petiole!(leaf, petiole, i, parameters["petiole_rachis_ratio_mean"], parameters["petiole_rachis_ratio_sd"], rng)
            for p in 2:nb_petiole_segments
                petiole = Node(petiole, NodeMTG("<", "Petiole", p, 5))
                compute_properties_petiole!(leaf, petiole, i, parameters["petiole_rachis_ratio_mean"], parameters["petiole_rachis_ratio_sd"], rng)
            end
            rachis = Node(petiole, NodeMTG("<", "Rachis", i, 5))
        end

        # add petiole, rachis, leaflets, ls
    end

    return plant
end
