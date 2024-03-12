"""
    read_param_file(file)

Reads a parameter file and returns the contents as an ordered dictionary.

# Arguments

- `file`: The path to the parameter file.

# Returns

An ordered dictionary containing the contents of the parameter file.

# Example

```julia
file = joinpath(dirname(dirname(pathof(VPalm))),"test","files","parameter_file.yml")
read_param_file(file)
```
"""
function read_param_file(file)
    YAML.load_file(file; dicttype=OrderedCollections.OrderedDict{String,Any})
end