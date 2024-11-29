"""
    read_parameters(file)

Reads a parameter file and returns the contents as an ordered dictionary.

# Arguments

- `file`: The path to the parameter file.

# Returns

An ordered dictionary containing the contents of the parameter file.

# Example

```julia
file = joinpath(dirname(dirname(pathof(VPalm))),"test","files","parameter_file.yml")
read_parameters(file)
```
"""
function read_parameters(file)
    p = YAML.load_file(file; dicttype=OrderedCollections.OrderedDict{String,Any})
    p["seed"] = p["seed"] |> Int
    p["nb_leaves_emitted"] = p["nb_leaves_emitted"] |> Int
    p["nb_leaves_mean"] = p["nb_leaves_mean"] |> Int
    p["nb_leaves_sd"] = p["nb_leaves_sd"] |> Int
    p["stem_growth_start"] = p["stem_growth_start"] |> Int
    p["nb_leaves_in_sheath"] = p["nb_leaves_in_sheath"] |> Int
    p["internode_rank_no_expansion"] = p["internode_rank_no_expansion"] |> Int
    p["nbInflorescences"] = p["nbInflorescences"] |> Int

    @assert p["nb_leaves_emitted"] > 0
    @assert p["nb_leaves_mean"] > 0
    return p
end


"""
    write_parameters(file, params)

Write the given parameters to a file using YAML format.

# Arguments
- `file`: The file path to write the parameters to.
- `params`: The parameters to be written.

# Example

```julia
file = joinpath(dirname(dirname(pathof(VPalm))),"test","files","parameter_file.yml")
params = read_parameters(file)
write_parameters(tempname(), params)
```
"""
function write_parameters(file, params)
    YAML.write_file(file, params)
end