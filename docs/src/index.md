```@meta
CurrentModule = VPalm
```

# VPalm

Documentation for [VPalm](https://github.com/PalmStudio/VPalm.jl).

Install VPalm using:

```julia
] add https://github.com/PalmStudio/VPalm.jl
```

Make a 3D palm using:

```julia
using VPalm # Import the package

# Import the example parameters:
parameters = read_parameters(joinpath(vpalm_test_files, "parameter_file.yml"), verbose=false)

# Build the 3D mockup:
palm = build_mockup(parameters)
```

Use `PlantGeom` and a `Makie` backend to visualize the palm:

```julia
using PlantGeom, CairoMakie
viz(palm)
```

## API

```@index
```

```@autodocs
Modules = [VPalm]
```
