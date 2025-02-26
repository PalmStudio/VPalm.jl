"""
    bend(
        type, width_bend, height_bend, init_torsion, x, y, z, mass_rachis, mass_leaflets_right, mass_leaflets_left,
        distance_application, elastic_modulus, shear_modulus, step, points, iterations;
        all_points=false,
        angle_max=deg2rad(21),
        force=true,
        verbose=true
    )

Compute the deformation of the rachis by applying both bending and torsion.

# Arguments
- `type`: Vector of section types (1: triangle bottom, 2: rectangle, 3: triangle top, 4: ellipse, 5: circle).
- `width_bend`: Vector of segment section widths (m).
- `height_bend`: Vector of segment section heights (m).
- `init_torsion`: Vector of initial torsion angles (degrees).
- `x`: Vector of x coordinates of the segments.
- `y`: Vector of y coordinates of the segments.
- `z`: Vector of z coordinates of the segments.
- `mass_rachis`: Vector of rachis segment masses (kg).
- `mass_leaflets_right`: Vector of leaflet masses carried by the segment, on the right side (kg).
- `mass_leaflets_left`: Vector of leaflet masses carried by the segment, on the left side (kg).
- `distance_application`: Vector of application distances for the left and right weights (m).
- `elastic_modulus`: Vector of elasticity moduli (bending, MPa).
- `shear_modulus`: Vector of shear moduli (torsion, MPa).
- `step`: Length of the segments that discretize the object (m).
- `points`: Number of points used in the grid discretizing the section.
- `iterations`: Number of iterations to compute the torsion and bending.
- `all_points=false`: return all points used in the computation (`true`), or only the input points corresponding to x, y and z coordinates (`false`, default).
- `angle_max=deg2rad(21)`: Maximum angle for testing the small displacement hypothesis (radians).
- `force=true`: Check if verify the small dispacements hypothesis and bounds the values to be at miximum `angle_max`
- `verbose=true`: Provide information during computation.

# Returns
Named tuple with geometrical fields describing the rachis bended and with torsion applied
- `x`: x coordinates of the points.
- `y`: y coordinates of the points.
- `z`: z coordinates of the points.
- `length`: length of the segments.
- `angle_xy`: angle between the xy-plan and the segment.
- `angle_xz`: angle between the xz-plan and the segment.
- `torsion`: torsion angle of the segment.
All these fields are vectors of the same length as the input vectors (i.e. number of segments).

# Details
The bending and torsion are applied to the sections of the rachis defined by 5 segments.

"""
function bend(type, width_bend, height_bend, init_torsion, x, y, z, mass_rachis, mass_leaflets_right, mass_leaflets_left,
    distance_application, elastic_modulus, shear_modulus, step, points, iterations;
    all_points=false,
    angle_max=deg2rad(21),
    force=true,
    verbose=true
)

    gravity = 9.8

    vec_rot_flex = zeros(3) # Originally a 3x1 matrix in R

    # Number of experimental points
    npoints_exp = length(x)  # Assuming x, y, z have the same length

    # Distances and angles of each segment P2P1
    vdist_p2p1, = xyz_to_dist_and_angles(x, y, z)

    dist_lineique = [0; cumsum(vdist_p2p1)] # For interpolation
    dist_totale = last(dist_lineique)

    # The distances of the segments cannot be zero. The origin point (0,0,0) cannot be in the data
    if any(vdist_p2p1 .== 0)
        error("Found distances between segments equal to 0.")
    end

    poids_lin_tige = mass_rachis ./ vdist_p2p1 ./ iterations
    v_poids_feuilles_d = mass_leaflets_right ./ vdist_p2p1 ./ iterations
    v_poids_feuilles_g = mass_leaflets_left ./ vdist_p2p1 ./ iterations
    v_poids_flexion = poids_lin_tige .+ v_poids_feuilles_d .+ v_poids_feuilles_g

    # Linear interpolations, segment length = step
    nlin = round(Int, dist_totale / step + 1)
    step = dist_totale / (nlin - 1)
    vec_dist = collect((0:(nlin-1)) .* step)
    vec_dist[end] = dist_lineique[end]
    # Note: we force vec_dist[end] to dist_lineique[end] to avoid any rounding error

    if length(elastic_modulus) != npoints_exp
        if length(elastic_modulus) == 1
            elastic_modulus = fill(elastic_modulus, npoints_exp)
        else
            error("`elastic_modulus` argument should be of length 1 or equal to `npoints_exp`")
        end
    end

    if length(shear_modulus) != npoints_exp
        if length(shear_modulus) == 1
            shear_modulus = fill(shear_modulus, npoints_exp)
        else
            error("`shear_modulus` argument should be of length 1 or equal to `npoints_exp`")
        end
    end

    vec_moe = linear_interpolation(dist_lineique, [elastic_modulus[1]; elastic_modulus])(vec_dist) .* 1e6
    vec_g = linear_interpolation(dist_lineique, [shear_modulus[1]; shear_modulus])(vec_dist) .* 1e6
    vangle_tor = deg2rad.(init_torsion)
    vec_agl_tor = linear_interpolation(dist_lineique, [vangle_tor[1]; vangle_tor])(vec_dist)
    vec_d_appli_poids_feuille = linear_interpolation(dist_lineique, [distance_application[1]; distance_application])(vec_dist)

    # Interpolation of coordinates in the origin frame
    # Identification of experimental points in the linear discretization
    vec_x, vec_y, vec_z, i_discret_pts_exp, vec_dist_p2p1, vec_angle_xy, vec_angle_xz = interp_points(x, y, z, step)

    val_epsilon = 1e-6
    if (vec_dist_p2p1[2] > (step + val_epsilon)) || (vec_dist_p2p1[2] < (step - val_epsilon))
        error("Point distance too narrow")
    end
    if (length(vec_x) != nlin)
        error("length(vec_x) != nlin")
    end

    # Increment of weight for the iterative calculation
    mat_dist_pts_exp = zeros(iterations, npoints_exp)

    som_cum_vec_agl_tor = copy(vec_agl_tor)  # geometric rotation and section

    for iter_poids in 1:iterations
        # Inertias and surfaces of the experimental points
        v_ig_flex = zeros(npoints_exp)
        v_ig_tor = zeros(npoints_exp)
        v_sr = zeros(npoints_exp)

        for iter in 1:npoints_exp
            ag_deg = rad2deg(som_cum_vec_agl_tor[i_discret_pts_exp[iter]])  # orientation section (degrees)
            inertia_flex_rot = inertia_flex_rota(width_bend[iter], height_bend[iter], ag_deg, type[iter], points) # Assuming this function is defined elsewhere
            v_ig_flex[iter] = inertia_flex_rot.ig_flex
            v_ig_tor[iter] = inertia_flex_rot.ig_tor
            v_sr[iter] = inertia_flex_rot.sr
        end

        # Linear interpolation of inertias
        vec_inertie_flex = linear_interpolation(dist_lineique, [v_ig_flex[1]; v_ig_flex])(vec_dist)
        vec_inertie_tor = linear_interpolation(dist_lineique, [v_ig_tor[1]; v_ig_tor])(vec_dist)

        # Write angles from the new coordinates
        # Distance and angles of each segment P2P1
        vec_dist_p2p1, vec_angle_xy, vec_angle_xz = xyz_to_dist_and_angles(vec_x, vec_y, vec_z) # Assuming this function is defined elsewhere

        vec_dist_p2p1[1] = 0
        vec_angle_xy[1] = vec_angle_xy[2]
        vec_angle_xz[1] = vec_angle_xz[2]

        # Flexion: linear bending forces and linear interpolation
        v_force = v_poids_flexion .* cos.(vec_angle_xy[i_discret_pts_exp]) .* gravity

        vec_force = linear_interpolation(dist_lineique, [v_force[1]; v_force])(vec_dist)

        # Shear forces and bending moments
        vec_shear = cumsum(vec_force[nlin:-1:1] .* step)
        vec_shear = vec_shear[nlin:-1:1]

        vec_moment = -cumsum(vec_shear[nlin:-1:1] .* step)
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
            force && (vec_angle_flexion[abs.(vec_angle_flexion).>angle_max] .= angle_max)
        end

        # Torsion
        v_m_tor = zeros(npoints_exp)

        for iter in 1:npoints_exp
            fdr = [0, 0, -v_poids_feuilles_d[iter] * vdist_p2p1[iter] * gravity]
            # Code with invariance by 'iterations'
            force_feuille_dr = rota_inverse_yz(fdr, vec_angle_xy[iter], vec_angle_xz[iter]) # Assuming this function is defined elsewhere

            fga = [0, 0, -v_poids_feuilles_g[iter] * vdist_p2p1[iter] * gravity]
            # Code with invariance by 'iterations'
            force_feuille_ga = rota_inverse_yz(fga, vec_angle_xy[iter], vec_angle_xz[iter]) # Assuming this function is defined elsewhere

            dist_point = vec_d_appli_poids_feuille[i_discret_pts_exp[iter]]
            angle_point = som_cum_vec_agl_tor[i_discret_pts_exp[iter]]

            # Hypothesis of contribution part right or left
            if angle_point > 0
                kd = 0
                kg = 1
            elseif angle_point < 0
                kd = 1
                kg = 0
            else # angle_point == 0
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
            force && (vec_angle_torsion[abs.(vec_angle_torsion).>angle_max] .= angle_max)
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
            vec_rot_inv = rota_inverse_yz(p2p1, vec_angle_xy[iter], vec_angle_xz[iter]) # Assuming this function is defined elsewhere

            # Flexion equivalent to a rotation around OY
            # Rotation around OY: The rotation is wrong for strong angles
            vec_rot_flex[1] = vec_rot_inv[1]
            vec_rot_flex[2] = vec_rot_inv[2]
            vec_rot_flex[3] = step * vec_angle_flexion[iter]

            # Torsion
            # Equivalent to a rotation around OX, initially the section is rotated but without torsion
            agl_tor_geom = som_cum_vec_agl_tor[iter] - vec_agl_tor[iter]

            cs = cos(agl_tor_geom)
            sn = sin(agl_tor_geom)

            mat_rot_x = [1 0 0; 0 cs -sn; 0 sn cs]

            vec_rot_tor = mat_rot_x * vec_rot_flex

            # Original base
            vec_rot = rota_yz(vec_rot_tor, vec_angle_xy[iter], vec_angle_xz[iter]) # Assuming this function is defined elsewhere

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
        XYZangles = xyz_to_dist_and_angles(vec_x, vec_y, vec_z)

        vec_x, vec_y, vec_z = dist_and_angles_to_xyz([0; fill(step, nlin - 1)], XYZangles.vangle_xy, XYZangles.vangle_xz) # Assuming this function is defined elsewhere

        # Calculation of the distances of the experimental points
        # Between before and after deformation
        for iter in 1:npoints_exp
            c1 = (x[iter] - vec_x[i_discret_pts_exp[iter]])^2
            c2 = (y[iter] - vec_y[i_discret_pts_exp[iter]])^2
            c3 = (z[iter] - vec_z[i_discret_pts_exp[iter]])^2

            mat_dist_pts_exp[iter_poids, iter] = sqrt(c1 + c2 + c3)
        end
    end

    if all_points
        i_discret_pts_exp = eachindex(vec_x)
    end

    pts_x = vec_x[i_discret_pts_exp]
    pts_y = vec_y[i_discret_pts_exp]
    pts_z = vec_z[i_discret_pts_exp]
    pts_agl_tor = rad2deg.(som_cum_vec_agl_tor[i_discret_pts_exp])

    pts_dist, pts_agl_xy, pts_agl_xz = xyz_to_dist_and_angles(pts_x, pts_y, pts_z) # Assuming this function is defined elsewhere

    return (x=pts_x, y=pts_y, z=pts_z, length=pts_dist, angle_xy=rad2deg.(pts_agl_xy), angle_xz=rad2deg.(pts_agl_xz), torsion=pts_agl_tor)
end
