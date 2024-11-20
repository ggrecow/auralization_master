function [ f ] = ita_diffraction_utd_illumination_fadeout( wedge, source_pos, receiver_pos, reflection_at_main_face )
%ita_diffraction_utd_illumination_fadeout Returns a scalar factor to fade out 
%last diffraction component between shadow and half way towards reflection boundary

assert( ~ita_diffraction_shadow_zone( wedge, source_pos, receiver_pos ) )
assert( ~ita_diffraction_reflection_zone( wedge, source_pos, receiver_pos, reflection_at_main_face ) )
%assert( ita_diffraction_utd_illumination_sign( wedge, source_pos, receiver_pos, reflection_at_main_face ) == -1 )


if reflection_at_main_face
    vw = wedge.main_face_normal;
else
    vw = wedge.opposite_face_normal;
end

vs = source_pos - wedge.location;
vr = receiver_pos - wedge.location;
distance_plane2source = dot( vs / norm( vs ) , vw / norm( vw ) );
distance_plane2receiver = dot( vr / norm( vr ) , vw / norm( vw ) );

if distance_plane2source == 0 || distance_plane2receiver == 0
    f = 1; % avoid NaNs in output, special case ...
    return
end

if reflection_at_main_face
    r = distance_plane2receiver / distance_plane2source;
    s = 1 - min( 1, abs( r * 1.1 ) ); % 0 ... 1 at 90 % through
    f = sin( s *  pi / 2 ) * sin( s *  pi / 2 ); % sine-square fading from 0 to 1
else
    r = distance_plane2source / distance_plane2receiver;
    s = 1 - max( 1, abs( r * 0.9 ) ); % 0 ... 1 at 90 % through
    f = cos( s *  pi / 2 ) * cos( s *  pi / 2 ); % cosine-square fading from 1 to 0
end

end
