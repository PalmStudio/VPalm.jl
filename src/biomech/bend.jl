function bend(data, elastic_modulus, shear_modulus, step=0.02, points=100, iterations=15, all_points=true, verbose=true)
    """
        bend(data, elastic_modulus, shear_modulus; step=0.02, points=100, iterations=15, all_points=true, verbose=true)

    Compute the deformation by applying both bending and torsion.

    # Arguments
    - `data::DataFrame`: Point data frame with each row being a point and each column being:
        - `x`: x coordinate of the segment
        - `y`: y coordinate of the segment
        - `z`: z coordinate of the segment
        - `type`: section type of the segment. 1: triangle (bottom-oriented); 2: rectangle; 3: triangle (top-oriented); 4: ellipsis; 5: circle.
        - `width`: segment section width (m)
        - `height`: segment section height (m)
        - `mass`: mass of the segment (kg)
        - `mass_right`: mass carried by the segment, on the right side (kg)
        - `mass_left`: mass carried by the segment, on the left side (kg)
        - `torsion`: initial torsion angle of the segment (degree)
        - `distance_application`: application distances for the left and right weights (m)
    - `elastic_modulus::Union{Float64, Vector{Float64}}`: Elasticity modulus (bending, MPa)
    - `shear_modulus::Union{Float64, Vector{Float64}}`: Shear modulus (torsion, MPa)
    - `step::Float64`: Length of the segments that discretize the object (m). Default is 0.02.
    - `points::Int`: Number of points used in the grid discretizing the section. Default is 100.
    - `iterations::Int`: Number of iterations to compute the torsion and bending. Default is 15.
    - `all_points::Bool`: Return all points or only the ones corresponding to the input data? Default is true.
    - `verbose::Bool`: Return information during process? Default is true.

    # Returns
    - `DataFrame`: A data frame with the following columns:
        - `x`: X coordinate (m)
        - `y`: Y coordinate (m)
        - `z`: Z coordinate (m)
        - `length`: length of the section (m)
        - `angle_xy`: angle of the segment with the XY plane (degree)
        - `angle_xz`: angle of the segment with the XZ plane (degree)
        - `torsion`: rotation angle of the segment (degree)

    # Examples
    ```julia
    filepath = system.file("extdata/6_EW01.22_17_kanan.txt", package="VPalm.jl")
    df = read_mat(filepath)
    # Un-bending the field measurements:
    df = unbend(df)

    # Adding the distance of application of the left and right weight:
    df.distance_application = distance_weight_sine(df.x)

    # (Re-)computing the deformation:
    df_bend = bend(df, elastic_modulus=2000, shear_modulus=400, step=0.02, points=100, iterations=15, verbose=true)
    ```
    """

    vecRotFlex = zeros(3, 1)

    NpointsExp = size(data, 1)

    vX = data[:, :x]
    vY = data[:, :y]
    vZ = data[:, :z]

    XYZangles = XYZ_to_angles(vX, vY, vZ)

    vDist_P2P1 = XYZangles.dist_P2P1
    vAngle_XY = XYZangles.vAngle_XY
    vAngle_XZ = XYZangles.vAngle_XZ

    distLineique = cumsum(vDist_P2P1)
    distTotale = distLineique[end]

    if any(vDist_P2P1 .== 0)
        error("Found distances between segments equal to 0.")
    end

    poidsTige = data[:, :mass]
    poidsFeuilD = data[:, :mass_right]
    poidsFeuilG = data[:, :mass_left]

    poidsLinTige = poidsTige ./ vDist_P2P1
    poidsLinFeuilD = poidsFeuilD ./ vDist_P2P1
    poidsLinFeuilG = poidsFeuilG ./ vDist_P2P1

    vPoidsFlexion = poidsLinTige + poidsLinFeuilD + poidsLinFeuilG
    vPoidsFeuillesD = poidsLinFeuilD
    vPoidsFeuillesG = poidsLinFeuilG

    Nlin = round(Int, distTotale / step + 1)
    step = distTotale / (Nlin - 1)
    vecDist = (0:(Nlin - 1)) * step

    if length(elastic_modulus) != size(data, 1)
        if length(elastic_modulus) == 1
            elastic_modulus = fill(elastic_modulus, size(data, 1))
        else
            error("elastic_modulus argument should be of length 1 or equal to `size(data, 1)`")
        end
    end

    if length(shear_modulus) != size(data, 1)
        if length(shear_modulus) == 1
            shear_modulus = fill(shear_modulus, size(data, 1))
        else
            error("shear_modulus argument should be of length 1 or equal to `size(data, 1)`")
        end
    end

    vMOE = elastic_modulus
    vecMOE = interpolate(c(0, distLineique), c(vMOE[1], vMOE), vecDist, method=:linear)

    vG = shear_modulus
    vecG = interpolate(c(0, distLineique), c(vG[1], vG), vecDist, method=:linear)

    vAngle_Tor = data[:, :torsion] * π / 180
    vecAglTor = interpolate(c(0, distLineique), c(vAngle_Tor[1], vAngle_Tor), vecDist, method=:linear)

    vDAppliPoidsFeuil = data[:, :distance_application]
    vecDAppliPoidsFeuil = interpolate(c(0, distLineique), c(vDAppliPoidsFeuil[1], vDAppliPoidsFeuil), vecDist, method=:linear)

    interp_list = InterpPoints(data, step)

    vecX = interp_list.vecX
    vecY = interp_list.vecY
    vecZ = interp_list.vecZ
    iDiscretPtsExp = interp_list.iDiscretPtsExp
    vecDist_P2P1 = interp_list.vecDist_P2P1
    vecAngle_XY = interp_list.vecAngle_XY
    vecAngle_XZ = interp_list.vecAngle_XZ

    valEpsilon = 1e-6
    if (vecDist_P2P1[2] > (step + valEpsilon)) || (vecDist_P2P1[2] < (step - valEpsilon))
        error("Point distance too narrow")
    end
    if length(vecX) != Nlin
        error("length(vecX) != Nlin")
    end

    vecMOE *= 1e6
    vecG *= 1e6

    matDistPtsExp = zeros(iterations, NpointsExp)

    vPoidsFlexion /= iterations
    vPoidsFeuillesD /= iterations
    vPoidsFeuillesG /= iterations

    SomCum_vecAglTor = vecAglTor

    for iterPoids in 1:iterations
        vIgFlex = vIgTor = vSr = zeros(NpointsExp)

        for iter in 1:NpointsExp
            b = data[iter, :width]
            h = data[iter, :height]
            sct = data[iter, :type]
            agDeg = SomCum_vecAglTor[iDiscretPtsExp[iter]] * 180 / π

            InertiaFlexRot = InertieFlexRota(b, h, agDeg, sct, points)
            IgFlex = InertiaFlexRot.IgFlex
            IgTor = InertiaFlexRot.IgTor
            Sr = InertiaFlexRot.Sr

            vIgFlex[iter] = IgFlex
            vIgTor[iter] = IgTor
            vSr[iter] = Sr
        end

        vecInertieFlex = interpolate(c(0, distLineique), c(vIgFlex[1], vIgFlex), vecDist, method=:linear)
        vecInertieTor = interpolate(c(0, distLineique), c(vIgTor[1], vIgTor), vecDist, method=:linear)

        XYZangles = XYZ_to_angles(vecX, vecY, vecZ)

        vecDist_P2P1 = XYZangles.dist_P2P1
        vecAngle_XY = XYZangles.vAngle_XY
        vecAngle_XZ = XYZangles.vAngle_XZ

        vecDist_P2P1[1] = 0
        vecAngle_XY[1] = vecAngle_XY[2]
        vecAngle_XZ[1] = vecAngle_XZ[2]

        pesanteur = 9.8

        vForce = vPoidsFlexion .* cos.(vecAngle_XY[iDiscretPtsExp]) .* pesanteur

        vecForce = interpolate(c(0, distLineique), c(vForce[1], vForce), vecDist, method=:linear)

        vecTranchant = cumsum(reverse(vecForce) .* step)
        vecTranchant = reverse(vecTranchant)

        vecMoment = -cumsum(reverse(vecTranchant) .* step)
        vecMoment = reverse(vecMoment)

        fct = vecMoment ./ (vecMOE .* vecInertieFlex)

        vecAngleFlexion = cumsum(reverse(fct) .* step)
        vecAngleFlexion = reverse(vecAngleFlexion)

        vecAngleFlexion = vecAngleFlexion[1] .- vecAngleFlexion

        AngleMax = 21 * π / 180

        if verbose && maximum(abs.(vecAngleFlexion)) > AngleMax
            println("Maximum bending angle (degree) = ", maximum(abs.(vecAngleFlexion)) * 180 / π)
            println("(!) Hypothesis of small displacements not verified for BENDING (!)")
        end

        vMTor = zeros(NpointsExp)

        for iter in 1:NpointsExp
            FDr = [0, 0, -vPoidsFeuillesD[iter] * vDist_P2P1[iter] * pesanteur]
            ForceFeuilDr = inverse_rotation_YZ(FDr, vecAngle_XY[iter], vecAngle_XZ[iter])

            FGa = [0, 0, -vPoidsFeuillesG[iter] * vDist_P2P1[iter] * pesanteur]
            ForceFeuilGa = inverse_rotation_YZ(FGa, vecAngle_XY[iter], vecAngle_XZ[iter])

            DistPoint = vecDAppliPoidsFeuil[iDiscretPtsExp[iter]]
            AnglePoint = SomCum_vecAglTor[iDiscretPtsExp[iter]]

            if AnglePoint > 0
                kD = 0
                kG = 1
            elseif AnglePoint < 0
                kD = 1
                kG = 0
            else
                kD = 0
                kG = 0
            end

            MD = DistPoint * kD * cos(AnglePoint) * ForceFeuilDr[3]
            MG = DistPoint * kG * cos(AnglePoint + π) * ForceFeuilGa[3]

            vMTor[iter] = MD + MG
        end

        vecMTor = interpolate(c(0, distLineique), c(vMTor[1], vMTor), vecDist, method=:linear)

        vecDerivAglTor = vecMTor ./ (vecG .* vecInertieTor)

        vecAngleTorsion = cumsum(vecDerivAglTor .* step)

        if maximum(abs.(vecAngleTorsion)) > AngleMax
            println("Maximum torsion angle (degree) = ", maximum(abs.(vecAngleTorsion)) * 180 / π)
            println("(!) Hypothesis of small displacements not verified for TORSION(!)")
        end

        SomCum_vecAglTor += vecAngleTorsion

        if verbose && iterPoids == iterations
            println(" ")
            println("Final torsion angle at the tip (degree) = ", SomCum_vecAglTor[end] * 180 / π)
        end

        neoVecX = neoVecY = neoVecZ = zeros(Nlin)

        for iter in 1:Nlin
            P2 = [vecX[iter], vecY[iter], vecZ[iter]]

            if iter == 1
                P1 = zeros(3)
            else
                P1 = [vecX[iter-1], vecY[iter-1], vecZ[iter-1]]
            end

            P2P1 = P2 - P1

            vecRotInv = inverse_rotation_YZ(P2P1, vecAngle_XY[iter], vecAngle_XZ[iter])

            vecRotFlex[1] = vecRotInv[1]
            vecRotFlex[2] = vecRotInv[2]
            vecRotFlex[3] = step * vecAngleFlexion[iter]

            aglTorGeom = SomCum_vecAglTor[iter] - vecAglTor[iter]

            cs = cos(aglTorGeom)
            sn = sin(aglTorGeom)

            matRotX = [1 0 0; 0 cs -sn; 0 sn cs]

            vecRotTor = matRotX * vecRotFlex

            vecRot = Rota_YZ(vecRotTor, vecAngle_XY[iter], vecAngle_XZ[iter])

            if iter == 1
                neoX = vecRot[1]
                neoY = vecRot[2]
                neoZ = vecRot[3]
            else
                neoX = neoVecX[iter - 1] + vecRot[1]
                neoY = neoVecY[iter - 1] + vecRot[2]
                neoZ = neoVecZ[iter - 1] + vecRot[3]
            end

            neoVecX[iter] = neoX
            neoVecY[iter] = neoY
            neoVecZ[iter] = neoZ
        end

        vecX = neoVecX
        vecY = neoVecY
        vecZ = neoVecZ

        XYZangles = XYZ_to_angles(vecX, vecY, vecZ)

        coords = angles_to_XYZ([0; fill(step, Nlin-1)], XYZangles.vAngle_XY, XYZangles.vAngle_XZ)

        vecX = coords.vecX
        vecY = coords.vecY
        vecZ = coords.vecZ

        for iter in 1:NpointsExp
            c1 = (vX[iter] - vecX[iDiscretPtsExp[iter]])^2
            c2 = (vY[iter] - vecY[iDiscretPtsExp[iter]])^2
            c3 = (vZ[iter] - vecZ[iDiscretPtsExp[iter]])^2

            matDistPtsExp[iterPoids, iter] = sqrt(c1 + c2 + c3)
        end
    end

    if all_points
        iDiscretPtsExp = 1:length(vecX)
    end

    PtsX = vecX[iDiscretPtsExp]
    PtsY = vecY[iDiscretPtsExp]
    PtsZ = vecZ[iDiscretPtsExp]

    Points = XYZ_to_angles(PtsX, PtsY, PtsZ)

    PtsDist = Points.dist_P2P1
    PtsAglXY = Points.vAngle_XY
    PtsAglXZ = Points.vAngle_XZ

    PtsAglXY = PtsAglXY * 180 / π
    PtsAglXZ = PtsAglXZ * 180 / π
    PtsAglTor = SomCum_vecAglTor[iDiscretPtsExp] * 180 / π

    return DataFrame(x=PtsX, y=PtsY, z=PtsZ, length=PtsDist, angle_xy=PtsAglXY, angle_xz=PtsAglXZ, torsion=PtsAglTor)
end