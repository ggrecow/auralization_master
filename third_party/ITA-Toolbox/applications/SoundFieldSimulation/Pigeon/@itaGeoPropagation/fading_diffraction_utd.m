function [ g_fading ] = fading_diffraction_utd( obj, anchor, source_pos, receiver_pos )
%fading_diffraction_utd Calculates the diffraction fading factor ( 0 .. +1 ) based on uniform
%theory of diffraction (illumination region, if no further reflection contribution is available)

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
if abs( cross(n1, n2) ) < itaInfiniteWedge.set_get_geo_eps
    w = itaSemiInfinitePlane( n1, loc, aperture_dir );
else
    w = itaInfiniteWedge( n1, n2, loc );
end

% Define section where reflection will occur

receiver_closer_main = ( dot( n1 - n2, receiver_pos - loc ) > 0 );
source_closer_main = ( dot( n1 - n2, source_pos - loc ) > 0 );

if receiver_closer_main && source_closer_main % but not in main face reflection zone ...
    reflection_at_main_face = true;
elseif ~receiver_closer_main && ~source_closer_main % but not in opposite face reflection zone ...
    reflection_at_main_face = false;
else
    % mixed situation
    nn = ( n1 + n2 ) / norm( n1 + n2 );
    sa = ( loc - source_pos ) / norm( loc - source_pos );
    ar = ( receiver_pos - loc ) / norm( receiver_pos - loc );
    nx = ( nn + sa ) / norm( nn + sa );
    if dot( nn, ar ) > dot( nn, nx )
        reflection_at_main_face = true;
    else
        reflection_at_main_face = false;
    end
end

if  ita_diffraction_shadow_zone( w, source_pos, receiver_pos )
    g_fading = 1; % if inside shadow zone, contribution fully requied
elseif ita_diffraction_reflection_zone( w, source_pos, receiver_pos, reflection_at_main_face )
    g_fading = 0; % if inside reflection zone, contribution must be removed
else
    g_fading = ita_diffraction_utd_illumination_fadeout( w, source_pos, receiver_pos, reflection_at_main_face ); % Fading window for opposite refl., closes roughly half-way
end

end
