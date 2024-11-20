function [ valid ] = validate_diffraction( obj, anchor, source_pos, receiver_pos )
%validate_diffraction Validates if diffraction calculation is applicable for given problem

if ~isfield( anchor, 'anchor_type' )
    error( 'The anchor argument does not contain a field "anchor_type"' )
end

valid = true;
    
% Assemble wedge
n1 = anchor.main_wedge_face_normal( 1:3 );
n2 = anchor.opposite_wedge_face_normal( 1:3 );
loc = anchor.vertex_start( 1:3 );
endPt = anchor.vertex_end( 1:3 );
len = norm( endPt - loc );
aperture_dir = ( endPt - loc ) / len;

% check if wedge is a screen
if abs( cross( n1, n2 ) ) < itaInfiniteWedge.set_get_geo_eps
    w = itaSemiInfinitePlane( n1, loc, aperture_dir );
else
    w = itaInfiniteWedge( n1, n2, loc );
end

if ~w.point_outside_wedge( source_pos ) || ...
   ~w.point_outside_wedge( receiver_pos )
    valid = false;
end

end
