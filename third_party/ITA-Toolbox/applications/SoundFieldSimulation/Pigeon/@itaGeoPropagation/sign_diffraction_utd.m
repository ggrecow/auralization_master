function [ g_sign ] = sign_diffraction_utd( obj, anchor, source_pos, receiver_pos )
%sign_diffraction_utd Calculates the diffraction sign (+1 or -1) based on uniform
%theory of diffraction (illumination region)

if ~isfield( anchor, 'anchor_type' )
    error( 'The anchor argument does not contain a field "anchor_type"' )
end

% Assemble wedge
n1 = anchor.main_wedge_face_normal( 1:3 );
n2 = anchor.opposite_wedge_face_normal( 1:3 );
loc = anchor.vertex_start( 1:3 );
endPt = anchor.vertex_end( 1:3 );
len = norm( endPt - loc );
aperture_dir = ( endPt - loc ) / len;

% check if wedge is a screen
if norm( cross(n1, n2) ) < itaInfiniteWedge.set_get_geo_eps
    w = itaSemiInfinitePlane( n1, loc, aperture_dir );
else
    w = itaInfiniteWedge( n1, n2, loc );
end

if ita_diffraction_shadow_zone( w, receiver_pos, source_pos )
    
    g_sign = +1;
    disp( '+1 Shadow region' )
    return
    
end

% Define section where reflection will occur

R = receiver_pos;
S = source_pos;

both_entities_on_same_side = false;

barrier = w.opening_angle / 2;
if ( w.angle_main_face( S ) <= barrier && w.angle_main_face( R ) <= barrier )
    reflection_at_main_face = true;
    both_entities_on_same_side = true;
elseif ( w.angle_main_face( S ) > barrier && w.angle_main_face( R ) > barrier )
    reflection_at_main_face = false;
    both_entities_on_same_side = true;
end

if both_entities_on_same_side

    if ita_diffraction_reflection_zone( w, R, S, reflection_at_main_face )
        g_sign = -1.0;
         fprintf('-1 Reflection zone (main face: %i)\n', reflection_at_main_face )
    else
        g_sign = +1.0;
        fprintf( '+1 Non-reflection zone (main face: %i)\n', reflection_at_main_face )
    end
    
    return
end

% mixed situation ... complicated

reflection_at_main_face = false;
flip_entities_symmetrical_solution = false;
if ( w.angle_main_face( S ) + w.angle_main_face( R ) ) / 2 < barrier
    reflection_at_main_face = true;
    if w.angle_main_face( S ) < w.angle_main_face( R )
        flip_entities_symmetrical_solution = true;
    end
else
    if w.angle_main_face( S ) > w.angle_main_face( R )
        flip_entities_symmetrical_solution = true;
    end
end

if flip_entities_symmetrical_solution
    g_sign = ita_diffraction_utd_illumination_sign(  w, R, S, reflection_at_main_face );
    fprintf( '%+1.0f illuminated, flipped (main face: %i)\n', g_sign, reflection_at_main_face )
else
    g_sign = ita_diffraction_utd_illumination_sign(  w, S, R, reflection_at_main_face );
    fprintf( '%+1.0f illuminated, non-flipped (main face: %i)\n', g_sign, reflection_at_main_face )
end

end
