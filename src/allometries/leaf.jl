

function compute_leaf_structure(rank)

end


"""
    leaf_insertion_angle(rank)

Compute the insertion angle of the leaf on the internode.

Note: The insertion angle is computed using a logistic function.

# Arguments

- `rank`: The rank of the leaf.
"""

function leaf_insertion_angle(rank)
    return logistic(rank, 90, 0.05, 40)
end