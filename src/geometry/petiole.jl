"""
    generate_petiole_mesh(params)

Create full petiole mesh based on parameters.

# Arguments


"""
function generate_petiole_mesh(petiole_length, width_base, height_base, width_cpoint, height_cpoint,
    insertion_angle, zenithal_cpoint_angle, azimuthal_angle, nb_sections=15)

    petiole_sections = properties_petiole_section(
        petiole_length, width_base, height_base, width_cpoint, height_cpoint,
        insertion_angle, zenithal_cpoint_angle, azimuthal_angle, nb_sections
    )
    # Create segments
    segments = []
    #! Type this vector

    for p in 1:petiole_sections
        # Create elliptical cross section
        petiole_segment = elliptical_cylinder(p.width / 2.0, p.height / 2.0, p.length) |> Meshes.Rotate(RotXY(p.zenithal_angle_local, p.azimuthal_angle_local))

        push!(segments, petiole_segment)
    end

    # Merge segments into final mesh
    final_mesh = merge(segments...)

    return final_mesh
end