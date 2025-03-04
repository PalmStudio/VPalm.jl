"""
    inertia_flex_rota(b, h, ag_deg, sct, n = 100)

Compute the inertia of bending and torsion, and the cross-section area.

# Arguments

- `b`: Dimension of the base.
- `h`: Dimension of the height.
- `ag_deg`: Section orientation angle (torsion, in degrees).
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
function inertia_flex_rota(b, h, ag_deg, sct, n=100)
    pas = min(b, h) / n
    nn = round(Int, h / pas)
    m = round(Int, b / pas)

    # Creation de la section
    section = zeros(nn, m)
    section = remplir(section, sct)

    # Construction des variables vectorisees
    iter_n = 1:size(section, 1)
    iter_m = 1:size(section, 2)

    mat_ind_ligne = repeat(iter_n, 1, m)
    mat_ind_colonne = repeat(iter_m', nn, 1)

    # Centre de gravite
    ng = sum(section .* mat_ind_ligne) / sum(section)
    mg = sum(section .* mat_ind_colonne) / sum(section)

    # Inerties et surface
    angle_radian = deg2rad(ag_deg)

    # Use Rotations.jl instead of manual rotation matrix
    rotation = RotZ(angle_radian)

    point_x = section .* ((mat_ind_colonne .- mg) .* pas)
    point_y = section .* ((mat_ind_ligne .- ng) .* pas)

    # Create 3D points (with z=0) for rotation
    points = [point_x[:] point_y[:] zeros(length(point_x[:]))]

    # Apply rotation to each point
    rotated_points = [rotation * [points[i, 1], points[i, 2], points[i, 3]] for i in 1:size(points, 1)]

    # Extract x and y coordinates from rotated points
    x = [p[1] for p in rotated_points]
    y = [p[2] for p in rotated_points]

    ds = pas^2
    ig_flex = sum(y .^ 2) * ds
    ig_tor = sum(x .^ 2 .+ y .^ 2) * ds
    sr = sum(section) * ds

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