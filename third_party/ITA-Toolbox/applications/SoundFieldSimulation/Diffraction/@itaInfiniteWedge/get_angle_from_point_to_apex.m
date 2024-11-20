function alpha_rad = get_angle_from_point_to_apex( obj, field_point, point_on_edge )
%Returns angle (radiant) between the ray from field point to apex point and the
%edge of the wedge.
%   output angle alpha: 0 <= alpha <= pi/2

if ~obj.point_outside_wedge( field_point )
    error( 'Field point must be outside wedge' );
end

if ~obj.point_on_edge( point_on_edge )
    warning( 'No point on edge found' )
    alpha_rad = [];
    return
end

dir_vec = ( point_on_edge - field_point ) / norm( point_on_edge - field_point );
alpha_rad = acos( dot( dir_vec, obj.edge_direction ) );

if alpha_rad > pi/2
    alpha_rad = pi - alpha_rad;
end

end
