function merged_outer_edge_diffraction = ita_propagation_merge_outer_edge_diffraction( a, b )
%ita_propagation_merge_anchors Merges two outer edge diffractions in a best-fit attempt

aperture_dir_a = a.vertex_end - a.vertex_start;
aperture_dir_b = b.vertex_end - b.vertex_start;

if a.main_wedge_face_id == b.main_wedge_face_id
    
    if dot( aperture_dir_a, aperture_dir_b ) > 0
        % Same direction, construct a centered substitute
    else
        % Opposite direction, use opposite b face
        merged_outer_edge_diffraction = a;
        merged_outer_edge_diffraction.main_wedge_face_normal = a.opposite_wedge_face_normal;
        merged_outer_edge_diffraction.main_wedge_face_id = a.opposite_wedge_face_id;
        merged_outer_edge_diffraction.opposite_wedge_face_normal = b.opposite_wedge_face_normal;
        merged_outer_edge_diffraction.opposite_wedge_face_id = b.opposite_wedge_face_id;
        merged_outer_edge_diffraction.vertex_start = mean( [ b.vertex_start, a.vertex_end ], 2 );
        merged_outer_edge_diffraction.vertex_end = mean( [ b.vertex_end, a.vertex_start ], 2 );
        merged_outer_edge_diffraction.interaction_point = mean( [ a.interaction_point, b.interaction_point ], 2 );
    end
    
elseif a.opposite_wedge_face_id == b.opposite_wedge_face_id
    
    if dot( aperture_dir_a, aperture_dir_b ) > 0
        % Same direction, construct a centered substitute
    else
        % Opposite direction, use opposite b face
        merged_outer_edge_diffraction = b;
        merged_outer_edge_diffraction.main_wedge_face_normal = b.opposite_wedge_face_normal;
        merged_outer_edge_diffraction.main_wedge_face_id = b.opposite_wedge_face_id;
        merged_outer_edge_diffraction.opposite_wedge_face_normal = a.opposite_wedge_face_normal;
        merged_outer_edge_diffraction.opposite_wedge_face_id = a.opposite_wedge_face_id;
        merged_outer_edge_diffraction.vertex_start = mean( [ a.vertex_start, b.vertex_end ], 2 );
        merged_outer_edge_diffraction.vertex_end = mean( [ a.vertex_end, b.vertex_start ], 2 );
        merged_outer_edge_diffraction.interaction_point = mean( [ b.interaction_point, a.interaction_point ], 2 );
    end

elseif a.main_wedge_face_id == b.opposite_wedge_face_id
    
    if dot( aperture_dir_a, aperture_dir_b ) > 0
        % Same direction, construct a centered substitute
        merged_outer_edge_diffraction = b;
        merged_outer_edge_diffraction.opposite_wedge_face_normal = a.opposite_wedge_face_normal;
        merged_outer_edge_diffraction.opposite_wedge_face_id = a.opposite_wedge_face_id;
        merged_outer_edge_diffraction.vertex_start = mean( [ a.vertex_start, b.vertex_start ], 2 );
        merged_outer_edge_diffraction.vertex_end = mean( [ a.vertex_end, b.vertex_end ], 2 );
        merged_outer_edge_diffraction.interaction_point = mean( [ a.interaction_point, b.interaction_point ], 2 );
    else
        % Opposite direction, use main b face
        merged_outer_edge_diffraction = b;
        merged_outer_edge_diffraction.opposite_wedge_face_normal = a.opposite_wedge_face_normal;
        merged_outer_edge_diffraction.opposite_wedge_face_id = a.opposite_wedge_face_id;
        merged_outer_edge_diffraction.vertex_start = mean( [ b.vertex_start, a.vertex_end ], 2 );
        merged_outer_edge_diffraction.vertex_end = mean( [ b.vertex_end, a.vertex_start ], 2 );
        merged_outer_edge_diffraction.interaction_point = mean( [ b.interaction_point, a.interaction_point ], 2 );
    end
    
elseif a.opposite_wedge_face_id == b.main_wedge_face_id
    
else
    % Not sharing geometric face
end

end
