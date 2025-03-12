"""
    unbend(df)

Removes torsion and bending of a bent beam by transforming it into a straight line,
while preserving its insertion angle (inclination angle of the first segment).

# Note
Mainly used to compute the input matrix for `bend()` from experimental points.

# Arguments
- `df`: DataFrame containing experimental data (see Details)

# Details
`df` must be a formatted DataFrame with each row representing a point and the following columns:
- `distance` (m): distance between previous point and this point (first value must be positive)
- `inclination` (degrees): insertion angle of the first point (only first value is used)

# Returns
The input DataFrame with modified (or enriched) columns:
- `x`: x coordinate
- `y`: y coordinate
- `z`: z coordinate

# Example
```julia
using CSV
filepath = joinpath(@__DIR__, "data", "6_EW01.22_17_kanan.txt")
df = CSV.read(filepath, DataFrame)
unbend(df)
```
"""
function unbend(df)
    # Distance between points cannot be 0
    zero_dist = findall(iszero, df.distance)
    if !isempty(zero_dist)
        df.distance[zero_dist] .= 1e-3 # (m)
    end
    # Cumulative distance of each segment
    x_distance = cumsum(df.distance)

    # Keep insertion angle of first segment
    agl_y = df.inclination[1] * π / 180
    agl_z = 0.0

    # Initialize columns if missing
    df.x = zeros(nrow(df))
    df.y = zeros(nrow(df))
    df.z = zeros(nrow(df))

    # Compute coordinates of points (unbent state)
    for iter in 1:nrow(df)
        op = [x_distance[iter], 0.0, 0.0]
        vec_rot = RotYZ(agl_y, agl_z) * op

        df.x[iter] = vec_rot[1]
        df.y[iter] = vec_rot[2]
        df.z[iter] = vec_rot[3]
    end

    return df
end