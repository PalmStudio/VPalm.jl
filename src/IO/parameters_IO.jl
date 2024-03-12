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
    YAML.load_file(file; dicttype=OrderedCollections.OrderedDict{String,Any})
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