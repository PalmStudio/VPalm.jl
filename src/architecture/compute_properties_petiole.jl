function compute_properties_petiole!(
    leaf_node,
    petiole_node,
    index,
    petiole_rachis_ratio_mean,
    petiole_rachis_ratio_sd,
    petiole_nb_sections,
    rng)

    leaf_node[:petioleLength] = mean_and_sd(petiole_rachis_ratio_mean, petiole_rachis_ratio_sd; rng=rng)
    petiole_node[:DeviationAngle] = normalDeviationDraw(5, rng)

    width_c_point = width_c_point(leaf_node[:rachis_length], c_point_width_intercept, c_point_width_slope)
    leaf_node[:c_point_angle] = c_point_angle(leaf_rank, c_point_decli_intercept, c_point_decli_slope, c_point_angle_SDP)

    petiole_segment_length = leaf_node[:petioleLength] / petiole_nb_sections
    petiole_node[:width] = petiole_width(width_c_point, petiole_section, petiole_nb_sections, leaf_base_width)
    petiole_node[:height] = petiole_height(height_c_point, petiole_section, petiole_nb_sections, leaf_base_height)
    petiole_node[:length] = petiole_segment_length

    return nothing
end






function width_c_point(rachis_length, c_point_width_intercept, c_point_width_slope)
    return linear(rachis_length, c_point_width_intercept, c_point_width_slope)
end


function c_point_angle(leaf_rank, c_point_decli_intercept, c_point_decli_slope, c_point_angle_sdp)
    angle = linear(leaf_rank, c_point_decli_intercept, c_point_decli_slope)
    angle += normal_deviation_draw(c_point_angle_sdp)
    angle = abs(angle)
    return leaf_rank < 3 ? 0.5 * angle : angle
end


function petiole_width(width_c_point, petiole_section, nb_petiole_sections, leaf_base_width)
    return leaf_base_width - (leaf_base_width - width_c_point) * (petiole_section / nb_petiole_sections)^0.17
end


function petiole_height(height_c_point, petiole_section, nb_petiole_sections, leaf_base_height)
    return leaf_base_height - (leaf_base_height - height_c_point) * (petiole_section / nb_petiole_sections)^0.5
end