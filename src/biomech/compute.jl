function XYZ_to_angles(vecX::Vector{Float64}, vecY::Vector{Float64}, vecZ::Vector{Float64})
    """
        XYZ_to_angles(vecX::Vector{Float64}, vecY::Vector{Float64}, vecZ::Vector{Float64}) -> NamedTuple

    Points to angles.

    # Arguments
    - `vecX::Vector{Float64}`: X coordinates (m)
    - `vecY::Vector{Float64}`: Y coordinates (m)
    - `vecZ::Vector{Float64}`: Z coordinates (m)

    # Returns
    - `NamedTuple`: A named tuple with the following fields:
        - `dist_P2P1`: Length of the segment
        - `vAngle_XY`: Angle between the segment and the XY plane (radian)
        - `vAngle_XZ`: Angle between the segment and the XZ plane (radian)

    # Examples
    ```julia
    filepath = joinpath(dirname(@__FILE__), "extdata", "6_EW01.22_17_kanan.txt")
    df = unbend(read_mat(filepath))
    XYZ_to_angles(df.x, df.y, df.z)
    ```
    """
    N = length(vecX)

    if length(vecY) != N
        error("Length of Y coordinates not equal to X coordinates")
    end
    if length(vecZ) != N
        error("Length of Z coordinates not equal to X coordinates")
    end

    dist_P2P1 = zeros(N)
    vAngle_XY = zeros(N)
    vAngle_XZ = zeros(N)

    for iter in 1:N
        P2 = [vecX[iter], vecY[iter], vecZ[iter]]

        if iter == 1
            P1 = [0.0, 0.0, 0.0]
        else
            P1 = [vecX[iter-1], vecY[iter-1], vecZ[iter-1]]
        end

        P2P1 = P2 .- P1

        if P2P1[1] != 0
            dist_P2P1[iter] = sqrt(P2P1[1]^2 + P2P1[2]^2 + P2P1[3]^2)
            vAngle_XY[iter] = atan(P2P1[3] / sqrt(P2P1[1]^2 + P2P1[2]^2))
            vAngle_XZ[iter] = atan(P2P1[2] / P2P1[1])
        end
    end

    return (dist_P2P1=dist_P2P1, vAngle_XY=vAngle_XY, vAngle_XZ=vAngle_XZ)
end

function angles_to_XYZ(dist_P2P1::Vector{Float64}, vAngle_XY::Vector{Float64}, vAngle_XZ::Vector{Float64})
    """
        angles_to_XYZ(dist_P2P1::Vector{Float64}, vAngle_XY::Vector{Float64}, vAngle_XZ::Vector{Float64}) -> NamedTuple

    Transform distances and angles into point coordinates.

    # Arguments
    - `dist_P2P1::Vector{Float64}`: Segment length (m)
    - `vAngle_XY::Vector{Float64}`: Angle between the segment and the XY plane (radian)
    - `vAngle_XZ::Vector{Float64}`: Angle between the segment and the XZ plane (radian)

    # Returns
    - `NamedTuple`: A named tuple with the following fields:
        - `vecX`: X coordinates
        - `vecY`: Y coordinates
        - `vecZ`: Z coordinates

    # Examples
    ```julia
    dist_P2P1 = [1.0, 2.0, 3.0]
    vAngle_XY = [0.1, 0.2, 0.3]
    vAngle_XZ = [0.1, 0.2, 0.3]
    angles_to_XYZ(dist_P2P1, vAngle_XY, vAngle_XZ)
    ```
    """
    N = length(dist_P2P1)

    if length(vAngle_XY) != N
        error("length of vAngle_XY != N")
    end
    if length(vAngle_XZ) != N
        error("length of vAngle_XZ != N")
    end

    vecX = zeros(N)
    vecY = zeros(N)
    vecZ = zeros(N)

    for iter in 1:N
        dZ = dist_P2P1[iter] * sin(vAngle_XY[iter])
        dist_XY = dist_P2P1[iter] * cos(vAngle_XY[iter])

        dX = dist_XY * cos(vAngle_XZ[iter])
        dY = dist_XY * sin(vAngle_XZ[iter])

        if iter == 1
            vecX[iter] = dX
            vecY[iter] = dY
            vecZ[iter] = dZ
        else
            vecX[iter] = vecX[iter-1] + dX
            vecY[iter] = vecY[iter-1] + dY
            vecZ[iter] = vecZ[iter-1] + dZ
        end
    end

    return (vecX=vecX, vecY=vecY, vecZ=vecZ)
end



function InterpPoints(matPoints::DataFrame, pas::Float64)
    """
        InterpPoints(matPoints::DataFrame, pas::Float64) -> NamedTuple

    Point interpolation.

    # Arguments
    - `matPoints::DataFrame`: Point matrix (from `unbend()`)
    - `pas::Float64`: Distance needed between interpolated points (m)

    # Returns
    - `NamedTuple`: A named tuple with the following fields:
        - `vecX`: Interpolated X coordinates
        - `vecY`: Interpolated Y coordinates
        - `vecZ`: Interpolated Z coordinates
        - `iDiscretPtsExp`: Indices of the experimental points in the interpolated points
        - `vecDist_P2P1`: Length of the segments
        - `vecAngle_XY`: Angle between the segments and the XY plane (radian)
        - `vecAngle_XZ`: Angle between the segments and the XZ plane (radian)
    """
    vX = matPoints[:, :x]
    vY = matPoints[:, :y]
    vZ = matPoints[:, :z]

    XYZ_Agl = XYZ_Vers_Agl(vX, vY, vZ)
    vDist_P2P1 = XYZ_Agl.dist_P2P1
    vAngle_XY = XYZ_Agl.vAngle_XY
    vAngle_XZ = XYZ_Agl.vAngle_XZ

    distLineique = cumsum(vDist_P2P1)
    distTotale = distLineique[end]

    if any(vDist_P2P1 .== 0)
        error("Found distances between segments equal to 0.")
    end

    Nlin = round(Int, distTotale / pas + 1)
    pas = distTotale / (Nlin - 1)

    NpointsExp = size(matPoints, 1)

    vecDist = [0; fill(pas, Nlin - 1)]
    distInterp = cumsum(vecDist)

    distLineique[end] += 1

    matXYZ = Matrix{Float64}(undef, 3, 0)

    for iter in 1:NpointsExp
        if iter == 1
            indPoints = distInterp .<= distLineique[iter]
        else
            indPoints = (distInterp .> distLineique[iter - 1]) .& (distInterp .<= distLineique[iter])
        end

        if !any(indPoints)
            error("No point found")
        end

        distPoints = vecDist[indPoints]
        OP = hcat(distPoints, zeros(length(distPoints)), zeros(length(distPoints)))'
        vecRot = Rota_YZ(OP, vAngle_XY[iter], vAngle_XZ[iter])

        if iter > 1
            vecPoints = cumsum(vecRot, dims=2) .+ matXYZ[:, end]
            matXYZ = hcat(matXYZ, vecPoints)
        else
            matXYZ = cumsum(vecRot, dims=2)
        end
    end

    vecX = matXYZ[1, :]
    vecY = matXYZ[2, :]
    vecZ = matXYZ[3, :]

    iDiscretPtsExp = zeros(Int, NpointsExp)

    for iter in 1:NpointsExp
        equad = sqrt.((vecX .- vX[iter]).^2 .+ (vecY .- vY[iter]).^2 .+ (vecZ .- vZ[iter]).^2)
        ind = argmin(equad)
        iDiscretPtsExp[iter] = ind
    end

    XYZ_Agl = XYZ_Vers_Agl(vecX, vecY, vecZ)
    vecDist_P2P1 = XYZ_Agl.dist_P2P1
    vecAngle_XY = XYZ_Agl.vAngle_XY
    vecAngle_XZ = XYZ_Agl.vAngle_XZ

    return (vecX=vecX, vecY=vecY, vecZ=vecZ, iDiscretPtsExp=iDiscretPtsExp, vecDist_P2P1=vecDist_P2P1, vecAngle_XY=vecAngle_XY, vecAngle_XZ=vecAngle_XZ)
end


function InertieFlexRota(b::Float64, h::Float64, agDeg::Float64, sct::Int, N::Int=100)
    """
        InertieFlexRota(b::Float64, h::Float64, agDeg::Float64, sct::Int, N::Int=100) -> NamedTuple

    Computes the inertia of bending and torsion, and the cross-section area.

    # Arguments
    - `b::Float64`: Dimension of the base
    - `h::Float64`: Dimension of the height
    - `agDeg::Float64`: Section orientation angle (torsion, in degrees)
    - `sct::Int`: Section type (see details)
    - `N::Int`: Number of discretizations (default to 100)

    # Details
    For the section type, possible values are:
    - `sct = 1`: triangle (bottom-oriented)
    - `sct = 2`: rectangle
    - `sct = 3`: triangle (top-oriented)
    - `sct = 4`: ellipsis
    - `sct = 5`: circle

    # Returns
    - `NamedTuple`: A named tuple with the following fields:
        - `IgFlex`: Bending inertia
        - `IgTor`: Torsion inertia
        - `Sr`: Cross-section surface
    """
    pas = min(b, h) / N
    n = round(Int, h / pas) + 1
    m = round(Int, b / pas) + 1

    section = zeros(n, m)
    section = Remplir(section, sct)

    iterN = 1:n
    iterM = 1:m

    matIndLigne = repeat(iterN', 1, length(iterM))
    matIndColonne = repeat(iterM, length(iterN), 1)

    ng = sum(section .* matIndLigne) / sum(section)
    mg = sum(section .* matIndColonne) / sum(section)

    angleRadian = agDeg * π / 180
    rotMatrix = [cos(angleRadian) -sin(angleRadian); sin(angleRadian) cos(angleRadian)]

    Point_x = section .* ((matIndColonne .- mg) .* pas)
    Point_y = section .* ((matIndLigne .- ng) .* pas)

    Point = [Point_x[:] Point_y[:]]'
    rotPoint = rotMatrix * Point

    x = rotPoint[1, :]
    y = rotPoint[2, :]

    dS = pas^2
    IgFlex = sum(y.^2) * dS
    IgTor = sum(x.^2 + y.^2) * dS
    Sr = sum(section) * dS

    return (IgFlex=IgFlex, IgTor=IgTor, Sr=Sr)
end


function inverse_rotation_YZ(OP::Matrix{Float64}, Agl_Y::Float64, Agl_Z::Float64)
    """
        inverse_rotation_YZ(OP::Matrix{Float64}, Agl_Y::Float64, Agl_Z::Float64) -> Matrix{Float64}

    Rotate around the Z axis and then the Y axis.

    # Arguments
    - `OP::Matrix{Float64}`: x, y, z coordinates
    - `Agl_Y::Float64`: Rotation angle around Y (radian)
    - `Agl_Z::Float64`: Rotation angle around Z (radian)

    # Returns
    - `Matrix{Float64}`: The rotated point

    # Examples
    ```julia
    inverse_rotation_YZ([0.0, 0.0, 0.0], 0.8517207, 0.0)
    ```
    """
    Agl_Y = -Agl_Y
    Agl_Z = -Agl_Z

    # Rotation around OZ
    cs = cos(Agl_Z)
    sn = sin(Agl_Z)

    matRotZ = [cs -sn 0; sn cs 0; 0 0 1]

    vecRotZ = matRotZ * OP

    # Rotation around OY
    cs = cos(Agl_Y)
    sn = sin(Agl_Y)

    matRotY = [cs 0 -sn; 0 1 0; sn 0 cs]

    return matRotY * vecRotZ
end