"""
    bend(type, width_bend, height_bend, init_torsion, x, y, z, mass, mass_right, mass_left,
         distance_application, elastic_modulus, shear_modulus, step, points, iterations)

Compute the deformation by applying both bending and torsion.

# Arguments
- `type`: Vector of section types (1: triangle bottom, 2: rectangle, 3: triangle top, 4: ellipse, 5: circle).
- `width_bend`: Vector of segment section widths (m).
- `height_bend`: Vector of segment section heights (m).
- `init_torsion`: Vector of initial torsion angles (degrees).
- `x`: Vector of x coordinates of the segments.
- `y`: Vector of y coordinates of the segments.
- `z`: Vector of z coordinates of the segments.
- `mass`: Vector of segment masses (kg).
- `mass_right`: Vector of masses carried by the segment, on the right side (kg).
- `mass_left`: Vector of masses carried by the segment, on the left side (kg).
- `distance_application`: Vector of application distances for the left and right weights (m).
- `elastic_modulus`: Vector of elasticity moduli (bending, MPa).
- `shear_modulus`: Vector of shear moduli (torsion, MPa).
- `step`: Length of the segments that discretize the object (m).
- `points`: Number of points used in the grid discretizing the section.
- `iterations`: Number of iterations to compute the torsion and bending.
- `all_points`: return all points used in the computation (`true`), or only the input points corresponding to x, y and z coordinates (`false`, default).
- `angle_max`: Maximum angle for testing the small displacement hypothesis (radians).
"""
function bend(type, width_bend, height_bend, init_torsion, x, y, z, mass, mass_right, mass_left,
    distance_application, elastic_modulus, shear_modulus, step, points, iterations;
    all_points=false,
    angle_max=deg2rad(21),  # 21 degrees is the limitGöttingen
    verbose=true
)

    vec_rot_flex = zeros(3) # Originally a 3x1 matrix in R

    # Number of experimental points
    npoints_exp = length(x)  # Assuming x, y, z have the same length

    # Distances and angles of each segment P2P1
    XYZangles = VPalm.xyz_vers_agl(x, y, z)

    vdist_p2p1 = XYZangles.dist_p2p1

    dist_lineique = [0; cumsum(vdist_p2p1)] # For interpolation 
    dist_totale = last(dist_lineique)

    # The distances of the segments cannot be zero
    # The origin point (0,0,0) cannot be in data
    if any(vdist_p2p1 .== 0)
        error("Found distances between segments equal to 0.")
    end

    # Linear weight of segments
    poids_tige = mass
    poids_feuille_d = mass_right
    poids_feuille_g = mass_left

    poids_lin_tige = poids_tige ./ vdist_p2p1
    poids_lin_feuille_d = poids_feuille_d ./ vdist_p2p1
    poids_lin_feuille_g = poids_feuille_g ./ vdist_p2p1

    v_poids_flexion = poids_lin_tige .+ poids_lin_feuille_d .+ poids_lin_feuille_g
    v_poids_feuilles_d = poids_lin_feuille_d
    v_poids_feuilles_g = poids_lin_feuille_g

    # Linear interpolations
    # Linear discretization
    # Segment length = step
    nlin = round(Int, dist_totale / step + 1)
    step = dist_totale / (nlin - 1)
    vec_dist = collect((0:(nlin-1)) .* step)

    if length(elastic_modulus) != npoints_exp
        if length(elastic_modulus) == 1
            elastic_modulus = fill(elastic_modulus, npoints_exp)
        else
            error("elastic_modulus argument should be of length 1 or equal to `npoints_exp`")
        end
    end

    if length(shear_modulus) != npoints_exp
        if length(shear_modulus) == 1
            shear_modulus = fill(shear_modulus, npoints_exp)
        else
            error("shear_modulus argument should be of length 1 or equal to `npoints_exp`")
        end
    end

    vec_moe = linear_interpolation(dist_lineique, [elastic_modulus[1]; elastic_modulus])(vec_dist)
    vec_g = linear_interpolation(dist_lineique, [shear_modulus[1]; shear_modulus])(vec_dist)
    vangle_tor = deg2rad.(init_torsion)
    vec_agl_tor = linear_interpolation(dist_lineique, [vangle_tor[1]; vangle_tor])(vec_dist)
    vec_d_appli_poids_feuille = linear_interpolation(dist_lineique, [distance_application[1]; distance_application])(vec_dist)

    # Interpolation of coordinates in the origin frame
    # Identification of experimental points in the linear discretization
    interp_list = VPalm.interp_points(x, y, z, step)

    vec_x = interp_list.vec_x
    vec_y = interp_list.vec_y
    vec_z = interp_list.vec_z
    i_discret_pts_exp = interp_list.i_discret_pts_exp
    vec_dist_p2p1 = interp_list.vec_dist_p2p1
    vec_angle_xy = interp_list.vec_angle_xy
    vec_angle_xz = interp_list.vec_angle_xz

    val_epsilon = 1e-6
    if (vec_dist_p2p1[2] > (step + val_epsilon)) || (vec_dist_p2p1[2] < (step - val_epsilon))
        error("Point distance too narrow")
    end
    if (length(vec_x) != nlin)
        error("length(vecX) != nlin")
    end

    # Increment of weight for the iterative calculation ===
    vec_moe = vec_moe .* 1e6
    vec_g = vec_g .* 1e6

    mat_dist_pts_exp = zeros(iterations, npoints_exp)

    v_poids_flexion = v_poids_flexion ./ iterations
    v_poids_feuilles_d = v_poids_feuilles_d ./ iterations
    v_poids_feuilles_g = v_poids_feuilles_g ./ iterations

    som_cum_vec_agl_tor = vec_agl_tor  # geometric rotation and section

    for iter_poids in 1:iterations

        # Inertias and surfaces of the experimental points
        v_ig_flex = zeros(npoints_exp)
        v_ig_tor = zeros(npoints_exp)
        v_sr = zeros(npoints_exp)

        for iter in 1:npoints_exp
            b = width_bend[iter]
            h = height_bend[iter]
            sct = type[iter]
            ag_deg = rad2deg(som_cum_vec_agl_tor[i_discret_pts_exp[iter]])  # orientation section (degrees)

            inertia_flex_rot = VPalm.inertia_flex_rota(b, h, ag_deg, sct, points) # Assuming this function is defined elsewhere
            ig_flex = inertia_flex_rot.ig_flex
            ig_tor = inertia_flex_rot.ig_tor
            sr = inertia_flex_rot.sr

            v_ig_flex[iter] = ig_flex
            v_ig_tor[iter] = ig_tor
            v_sr[iter] = sr
        end

        # Linear interpolation of inertias
        vec_inertie_flex = linear_interpolation(dist_lineique, [v_ig_flex[1]; v_ig_flex])(vec_dist)
        vec_inertie_tor = linear_interpolation(dist_lineique, [v_ig_tor[1]; v_ig_tor])(vec_dist)

        # Write angles from the new coordinates
        # Distance and angles of each segment P2P1
        XYZangles = VPalm.xyz_vers_agl(vec_x, vec_y, vec_z) # Assuming this function is defined elsewhere

        vec_dist_p2p1 = XYZangles.dist_p2p1
        vec_angle_xy = XYZangles.vangle_xy
        vec_angle_xz = XYZangles.vangle_xz

        vec_dist_p2p1[1] = 0
        vec_angle_xy[1] = vec_angle_xy[2]
        vec_angle_xz[1] = vec_angle_xz[2]

        # Flexion
        # Linear bending forces
        # and Linear interpolation
        pesanteur = 9.8

        # Linear force
        # Code with invariance by 'iterations'
        v_force = v_poids_flexion .* cos.(vec_angle_xy[i_discret_pts_exp]) .* pesanteur

        vec_force = linear_interpolation(dist_lineique, [v_force[1]; v_force])(vec_dist)

        # Shear forces and bending moments
        vec_tranchant = cumsum(vec_force[nlin:-1:1] .* step)
        vec_tranchant = vec_tranchant[nlin:-1:1]

        vec_moment = -cumsum(vec_tranchant[nlin:-1:1] .* step)
        vec_moment = vec_moment[nlin:-1:1]

        # Classic calculation of the deflection (distance delta)
        fct = vec_moment ./ (vec_moe .* vec_inertie_flex)

        vec_angle_flexion = cumsum(fct[nlin:-1:1] .* step)
        vec_angle_flexion = vec_angle_flexion[nlin:-1:1]

        # Embedded condition (derivative 1 = 0)
        vec_angle_flexion = vec_angle_flexion[1] .- vec_angle_flexion

        # Test of the small displacement hypothesis
        if verbose && maximum(abs.(vec_angle_flexion)) > angle_max
            @warn string("Maximum bending angle: ", rad2deg(maximum(abs.(vec_angle_flexion))), "°. Hypothesis of small displacements not verified for bending.")
        end

        # Torsion
        v_m_tor = zeros(npoints_exp)

        for iter in 1:npoints_exp
            fdr = [0, 0, -v_poids_feuilles_d[iter] * vdist_p2p1[iter] * pesanteur]
            # Code with invariance by 'iterations'
            force_feuille_dr = VPalm.rota_inverse_yz(fdr, vec_angle_xy[iter], vec_angle_xz[iter]) # Assuming this function is defined elsewhere

            fga = [0, 0, -v_poids_feuilles_g[iter] * vdist_p2p1[iter] * pesanteur]
            # Code with invariance by 'iterations'
            force_feuille_ga = VPalm.rota_inverse_yz(fga, vec_angle_xy[iter], vec_angle_xz[iter]) # Assuming this function is defined elsewhere

            dist_point = vec_d_appli_poids_feuille[i_discret_pts_exp[iter]]
            angle_point = som_cum_vec_agl_tor[i_discret_pts_exp[iter]]

            # Hypothesis of contribution part D or G
            if angle_point > 0
                kd = 0
                kg = 1
            elseif angle_point < 0
                kd = 1
                kg = 0
            elseif angle_point == 0
                kd = 0
                kg = 0
            end

            md = dist_point * kd * cos(angle_point) * force_feuille_dr[3]
            mg = dist_point * kg * cos(angle_point + pi) * force_feuille_ga[3]

            v_m_tor[iter] = md + mg
        end

        vec_m_tor = linear_interpolation(dist_lineique, [v_m_tor[1]; v_m_tor])(vec_dist)

        vec_deriv_agl_tor = vec_m_tor ./ (vec_g .* vec_inertie_tor)

        vec_angle_torsion = cumsum(vec_deriv_agl_tor .* step)  # integration along the stem

        if verbose && maximum(abs.(vec_angle_torsion)) > angle_max
            @warn string("Maximum torsion angle: ", rad2deg(maximum(abs.(vec_angle_torsion))), "°. Hypothesis of small displacements not verified for torsion.")
        end

        som_cum_vec_agl_tor = som_cum_vec_agl_tor .+ vec_angle_torsion  # cumulative by weight increment

        if verbose && iter_poids == iterations
            @info string("Final torsion angle at the tip: ", rad2deg(som_cum_vec_agl_tor[length(som_cum_vec_agl_tor)]), "°")
        end

        # New coordinates of the points
        neo_vec_x = zeros(nlin)
        neo_vec_y = zeros(nlin)
        neo_vec_z = zeros(nlin)

        for iter in 1:nlin
            # Origin P1
            p2 = [vec_x[iter], vec_y[iter], vec_z[iter]]

            if iter == 1
                p1 = zeros(3)
            else
                p1 = [vec_x[iter-1], vec_y[iter-1], vec_z[iter-1]]
            end

            p2p1 = p2 .- p1

            # Change of basis
            # Segment becomes collinear to the OX axis
            vec_rot_inv = VPalm.rota_inverse_yz(p2p1, vec_angle_xy[iter], vec_angle_xz[iter]) # Assuming this function is defined elsewhere

            # Flexion
            # Equivalent to a rotation around OY
            # Rotation around OY
            # The rotation is wrong for strong angles
            # Code replaced by:
            vec_rot_flex[1] = vec_rot_inv[1]
            vec_rot_flex[2] = vec_rot_inv[2]
            vec_rot_flex[3] = step * vec_angle_flexion[iter]

            # Torsion
            # Equivalent to a rotation around OX
            # Initially the section is rotated but without torsion
            agl_tor_geom = som_cum_vec_agl_tor[iter] - vec_agl_tor[iter]

            cs = cos(agl_tor_geom)
            sn = sin(agl_tor_geom)

            mat_rot_x = [1 0 0; 0 cs -sn; 0 sn cs]

            vec_rot_tor = mat_rot_x * vec_rot_flex

            # Original base
            vec_rot = VPalm.rota_yz(vec_rot_tor, vec_angle_xy[iter], vec_angle_xz[iter]) # Assuming this function is defined elsewhere

            # Point in the origin base
            if iter == 1
                neo_x = vec_rot[1]
                neo_y = vec_rot[2]
                neo_z = vec_rot[3]
            else
                neo_x = neo_vec_x[iter-1] + vec_rot[1]
                neo_y = neo_vec_y[iter-1] + vec_rot[2]
                neo_z = neo_vec_z[iter-1] + vec_rot[3]
            end

            # Re-writing the points
            neo_vec_x[iter] = neo_x
            neo_vec_y[iter] = neo_y
            neo_vec_z[iter] = neo_z
        end

        # Update variables
        vec_x = neo_vec_x
        vec_y = neo_vec_y
        vec_z = neo_vec_z

        # Conservation of distances
        # step = distance between points
        XYZangles = VPalm.xyz_vers_agl(vec_x, vec_y, vec_z) # Assuming this function is defined elsewhere

        coords = VPalm.agl_vers_xyz([0; fill(step, nlin - 1)], XYZangles.vangle_xy, XYZangles.vangle_xz) # Assuming this function is defined elsewhere

        vec_x = coords.vec_x
        vec_y = coords.vec_y
        vec_z = coords.vec_z

        # Calculation of the distances of the experimental points
        # Between before and after deformation
        for iter in 1:npoints_exp
            c1 = (x[iter] - vec_x[i_discret_pts_exp[iter]])^2
            c2 = (y[iter] - vec_y[i_discret_pts_exp[iter]])^2
            c3 = (z[iter] - vec_z[i_discret_pts_exp[iter]])^2

            mat_dist_pts_exp[iter_poids, iter] = sqrt(c1 + c2 + c3)
        end
    end  # iterPoids

    if all_points
        i_discret_pts_exp = eachindex(vec_x)
    end
    pts_x = vec_x[i_discret_pts_exp]
    pts_y = vec_y[i_discret_pts_exp]
    pts_z = vec_z[i_discret_pts_exp]

    Points = VPalm.xyz_vers_agl(pts_x, pts_y, pts_z) # Assuming this function is defined elsewhere

    pts_dist = Points.dist_p2p1
    pts_agl_xy = Points.vangle_xy
    pts_agl_xz = Points.vangle_xz

    pts_agl_xy = rad2deg.(pts_agl_xy)
    pts_agl_xz = rad2deg.(pts_agl_xz)
    pts_agl_tor = rad2deg.(som_cum_vec_agl_tor[i_discret_pts_exp])

    return (x=pts_x, y=pts_y, z=pts_z, length=pts_dist, angle_xy=pts_agl_xy, angle_xz=pts_agl_xz, torsion=pts_agl_tor)
end
