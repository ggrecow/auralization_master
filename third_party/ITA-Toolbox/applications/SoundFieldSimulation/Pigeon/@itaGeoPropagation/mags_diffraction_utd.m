function [ mags_linear, valid ] = mags_diffraction_utd( obj, anchor, source_pos, receiver_pos, f )
%mags_diffraction_utd Calculates the diffraction mgnitudes based on uniform
%theory of diffraction (with Kawai approximation). 

valid = false;
mags_linear = zeros( numel( f ), 1 );

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

% Validate
apex_point = w.get_aperture_point_far_field( source_pos, receiver_pos );
if all( apex_point == source_pos  ) || ...
   all ( apex_point == source_pos )
    warning( 'Skipping a path segment with a double anchor point' )
    return
elseif any( isnan( apex_point ) )
    warning( 'Skipping a path segment with an invalid aperture point' )
    return
end

apx = w.apex_point_approx( source_pos, receiver_pos );
[ ~, D, A ] = ita_diffraction_utd( w, source_pos, receiver_pos, f, obj.c, apx );

mags_linear = abs( A .* D );
valid = true;

end
