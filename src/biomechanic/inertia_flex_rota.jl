"""
    inertia_flex_rota(b, h, orientation_angle, sct, n = 100)

Compute the inertia of bending and torsion, and the cross-section area.

# Arguments

- `b`: Dimension of the base.
- `h`: Dimension of the height.
- `orientation_angle`: Section orientation angle (torsion, in radians).
- `sct`: Section type (see details).
- `n`: Number of discretizations (default to 100).

# Details

For the section type, possible values are:
- `sct = 1`: triangle (bottom-oriented)
- `sct = 2`: rectangle
- `sct = 3`: triangle (top-oriented)
- `sct = 4`: ellipse
- `sct = 5`: circle

# Returns

- A NamedTuple with fields:
  - `ig_flex`: Bending inertia.
  - `ig_tor`: Torsion inertia.
  - `sr`: Cross-section surface.
"""
function inertia_flex_rota(b, h, orientation_angle, sct, n=100)
    pas = min(b, h) / n
    nn = round(Int, h / pas)
    m = round(Int, b / pas)

    # Creation de la section
    section = zeros(nn, m)
    section = VPalm.remplir(section, sct)

    # Centre de gravite
    total_mass = sum(section)
    center_y = sum(section .* (1:nn)) / total_mass
    center_x = sum(section .* (1:m)') / total_mass

    # Create Meshes.Point collection
    p_type = typeof(Meshes.Point(0.0, 0.0, 0.0))
    points = Vector{p_type}()

    for i in 1:nn
        for j in 1:m
            if section[i, j] > 0
                # Calculate position relative to center of mass
                x = (j - center_x) * pas
                y = (i - center_y) * pas
                z = 0.0
                push!(points, Meshes.Point(x, y, z))
            end
        end
    end

    # Apply rotation to each point
    rotation = RotZ(orientation_angle)
    rotated_points = Meshes.Rotate(rotation)(points)

    # Calculate inertias
    ds = pas^2
    ig_flex = sum(Meshes.coords(p).y^2 for p in rotated_points) * ds
    ig_tor = sum(Meshes.coords(p).x^2 + Meshes.coords(p).y^2 for p in rotated_points) * ds
    sr = length(points) * ds

    return (ig_flex=ig_flex, ig_tor=ig_tor, sr=sr)
end

"""
    remplir(section, sct)

Fill in the matrix according to the section shape.

# Arguments
- `section`: Section matrix.
- `sct`: Section type (see details).

# Returns
- The filled matrix.
"""
function remplir(section, sct)
    # sct = 1 : triangle bas
    if sct == 1
        b13 = [1 1; size(section, 2)/2 1] \ [1; size(section, 1)]
        b23 = [size(section, 2) 1; size(section, 2)/2 1] \ [1; size(section, 1)]

        iter_n = 1:size(section, 1)
        iter_m = 1:size(section, 2)

        mat_ind_ligne = repeat(iter_n, 1, size(section, 2))
        mat_ind_colonne = repeat(iter_m', size(section, 1), 1)

        n13 = mat_ind_colonne * b13[1] .+ b13[2]
        n23 = mat_ind_colonne * b23[1] .+ b23[2]

        section = (mat_ind_ligne .<= n13) .& (mat_ind_ligne .<= n23)
    end

    # sct = 2 : rectangle
    if sct == 2
        section = ones(size(section))
    end

    # sct = 3 : triangle haut
    if sct == 3
        b13 = [1 1; size(section, 2)/2 1] \ [size(section, 1); 1]
        b23 = [size(section, 2) 1; size(section, 2)/2 1] \ [size(section, 1); 1]

        iter_n = 1:size(section, 1)
        iter_m = 1:size(section, 2)

        mat_ind_ligne = repeat(iter_n, 1, size(section, 2))
        mat_ind_colonne = repeat(iter_m', size(section, 1), 1)

        n13 = mat_ind_colonne * b13[1] .+ b13[2]
        n23 = mat_ind_colonne * b23[1] .+ b23[2]

        section = (mat_ind_ligne .>= n13) .& (mat_ind_ligne .>= n23)
    end

    # sct = 4 : ellipse
    if sct == 4
        a = maximum(size(section)) / 2
        b = minimum(size(section)) / 2

        c = sqrt(a^2 - b^2)

        iter_n = 1:size(section, 1)
        iter_m = 1:size(section, 2)

        mat_ind_ligne = repeat(iter_n, 1, size(section, 2))
        mat_ind_colonne = repeat(iter_m', size(section, 1), 1)

        if size(section, 1) >= size(section, 2)
            mf = size(section, 2) / 2

            nf1 = (a - c)
            nf2 = 2 * c + (a - c)

            dist1 = sqrt.((mat_ind_ligne .- nf1) .^ 2 .+ (mat_ind_colonne .- mf) .^ 2)
            dist2 = sqrt.((mat_ind_ligne .- nf2) .^ 2 .+ (mat_ind_colonne .- mf) .^ 2)

            section = ((dist1 .+ dist2) .<= (2 * a))
        end

        if size(section, 1) < size(section, 2)
            nf = size(section, 1) / 2

            mf1 = (a - c)
            mf2 = 2 * c + (a - c)

            dist1 = sqrt.((mat_ind_ligne .- nf) .^ 2 .+ (mat_ind_colonne .- mf1) .^ 2)
            dist2 = sqrt.((mat_ind_ligne .- nf) .^ 2 .+ (mat_ind_colonne .- mf2) .^ 2)

            section = ((dist1 .+ dist2) .<= (2 * a))
        end
    end

    # sct = 5 : cercle
    if sct == 5
        rayon = minimum(size(section)) / 2

        n0 = size(section, 1) / 2
        m0 = size(section, 2) / 2

        iter_n = 1:size(section, 1)
        iter_m = 1:size(section, 2)

        mat_ind_ligne = repeat(iter_n, 1, size(section, 2))
        mat_ind_colonne = repeat(iter_m', size(section, 1), 1)

        dist = sqrt.((mat_ind_ligne .- n0) .^ 2 .+ (mat_ind_colonne .- m0) .^ 2)
        section = dist .<= rayon
    end
    return section
end