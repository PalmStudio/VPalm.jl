function unbend(df::DataFrame)
    """
        unbend(df::DataFrame) -> DataFrame

    Remove torsion and bending of a bent beam, i.e. transform it into a straight line, while keeping its insertion angle (inclination angle of the first segment).

    # Arguments
    - `df::DataFrame`: Experimental data (usually read using `read_mat()`, see details)

    # Details
    `df` should be a formatted DataFrame, with each row being a point and each column being:
    - `distance`: distance between the previous point and this point (first value should be positive) (m)
    - `inclination`: insertion angle of the first point (only first value is used) (degree)

    # Returns
    - `DataFrame`: The input DataFrame with modified (or enriched with):
        - `x`: x coordinate
        - `y`: y coordinate
        - `z`: z coordinate

    # Examples
    ```julia
    filepath = joinpath(dirname(@__FILE__), "extdata", "6_EW01.22_17_kanan.txt")
    df = read_mat(filepath)
    unbend(df)
    ```
    """
    iD0 = findall(x -> x == 0, df.distance)
    if !isempty(iD0)
        df.distance[iD0] .= 1e-3 # (m)
    end

    Distance = df.distance

    XDistance = cumsum(Distance)

    Agl_Y = df.inclination[1] * π / 180
    Agl_Z = 0.0

    df.x = zeros(size(df, 1))
    df.y = zeros(size(df, 1))
    df.z = zeros(size(df, 1))

    for iter in 1:size(df, 1)
        OP = [XDistance[iter], 0.0, 0.0]
        vecRot = Rota_YZ(OP, Agl_Y, Agl_Z)

        df.x[iter] = vecRot[1]
        df.y[iter] = vecRot[2]
        df.z[iter] = vecRot[3]
    end

    return df
end