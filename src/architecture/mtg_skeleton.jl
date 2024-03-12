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

- `nb_leaves_emitted`: The number of leaves to emit.

# Examples

```julia
mtg_skeleton(3)
```
"""
function mtg_skeleton(nb_leaves_emitted)
    plant = Node(NodeMTG("/", "Plant", 1, 1))
    #roots = Node(plant, NodeMTG("+", "RootSystem", 1, 2))
    stem = Node(plant, NodeMTG("+", "Stem", 1, 2))

    phytomer = Node(stem, NodeMTG("/", "Phytomer", 1, 3))
    internode = Node(phytomer, NodeMTG("/", "Internode", 1, 4))
    Node(internode, NodeMTG("+", "Leaf", 1, 4))

    for i in 2:nb_leaves_emitted
        phytomer = Node(phytomer, NodeMTG("<", "Phytomer", i, 3))
        internode = Node(phytomer, NodeMTG("/", "Internode", i, 4))
        Node(internode, NodeMTG("+", "Leaf", i, 4))
    end

    return plant
end