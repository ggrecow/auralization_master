function [ diffr_comp_sign ] = ita_diffraction_utd_illumination_sign( wedge, source_pos, receiver_pos, reflection_at_main_face )
%ITA_DIFFRACTION_UTD_ILLUMINATION_SIGN Returns -1 if source is closer to
%shadow boundary and 1 if source is closer to reflection boundary 

if ita_diffraction_shadow_zone( wedge, source_pos, receiver_pos )
    warning 'Sign requested inside shadow zone, always positive here.'
    diffr_comp_sign = 1;
    return
end
if ita_diffraction_reflection_zone( wedge, source_pos, receiver_pos, reflection_at_main_face )
    warning 'Sign requested inside reflection zone, always negative here.'
    diffr_comp_sign = -1;
    return
end

vs = source_pos - wedge.location;
vr = receiver_pos - wedge.location;

if reflection_at_main_face
    vw = wedge.main_face_normal;
else
    vw = wedge.opposite_face_normal;
end

distance_plane2source = dot( vs / norm( vs ) , vw / norm( vw ) );
distance_plane2receiver = dot( vr / norm( vr ) , vw / norm( vw ) );

if distance_plane2source > 0
    diffr_comp_sign = 1; % closer to reflection boundary
else
    diffr_comp_sign = -1; % closer to shadow boundary
end

end
