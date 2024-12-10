
"""
    read_mat(path::String) -> DataFrame

Read the matrix from field measurements.

# Arguments
- `path::String`: Path to the file.

# Returns
- `DataFrame`: A formatted DataFrame, with each row being a point and with columns being:
    - `distance`: distance between the previous point and this point (first value should be positive) (m)
    - `type`: section type. 1: triangle (bottom-oriented); 2: rectangle; 3: triangle (top-oriented); 4: ellipsis; 5: circle.
    - `width`: section width (m)
    - `height`: section height (m)
    - `inclination`: insertion angle of the first point (only first value is used) (degree)
    - `torsion`: torsion angle (degree)
    - `x`: x coordinate
    - `y`: y coordinate
    - `z`: z coordinate
    - `mass`: mass of the section at the point (kg)
    - `mass_right`: mass carried by the object, on the right side (kg)
    - `mass_left`: mass carried by the object, on the left side (kg)

# Examples
```julia
filepath = joinpath(dirname(@__FILE__), "extdata", "6_EW01.22_17_kanan.txt")
df = read_mat(filepath)
```
"""

function read_mat(path::String)::DataFrame
    df = CSV.read(path, DataFrame)
    if ncol(df) != 12
        error("File not in the right format, expected 12 columns, found $(ncol(df)).")
    end
    df = DataFrame(transpose(Matrix(df)))
    rename!(df, [:distance, :type, :width, :height, :inclination, :torsion, :x, :y, :z, :mass, :mass_right, :mass_left])
    return df
end