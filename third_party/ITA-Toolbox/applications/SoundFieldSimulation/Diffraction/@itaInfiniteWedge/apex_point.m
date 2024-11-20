function ap = apex_point( obj, source_pos, receiver_pos )
% GET_APERTURE_POINT_FAR_FIELD Returns aperture point on wedge (closest point on wedge
% between source and receiver if both are in the far field)
    
assert( numel( source_pos ) == 3 )
assert( numel( receiver_pos ) == 3 )


%% Calculations

source_receiver_vec = receiver_pos - source_pos;
edge_dir = obj.edge_direction / norm( obj.edge_direction );
    
if norm( source_receiver_vec ) ~= 0

    % Based on a line-plane intersection
    source_receiver_dir = source_receiver_vec / norm( source_receiver_vec );

    % Auxilary plane spanned by source_receiver_dir and aux_plane_dir (closest
    % line between edge vector and source-receiver-vector
    aux_plane_vec = cross( source_receiver_dir, edge_dir );

    if norm( aux_plane_vec ) ~= 0

        % Directions are not parallel, a closest point must exist
        aux_plane_dir = aux_plane_vec / norm( aux_plane_vec );
        aux_plane_normal = cross( source_receiver_dir, aux_plane_dir );

        % Determine intersection of line (aperture) and auxiliary plane
        lambda_divisor = dot( aux_plane_normal, edge_dir );
        assert( lambda_divisor ~= 0 )
        d = dot( aux_plane_normal, obj.location ); % Distance to origin
        lambda = ( d - dot( aux_plane_normal, source_pos ) );% / lambda_divisor; % .. or receiver_pos
        ap = obj.location + lambda * edge_dir;

    else

        % Directions are parallel, project middle point between source & target
        % onto aperture

        % Project middle point onto aperture
        lambda = dot( ( source_pos + source_receiver_vec ./ 2 ) - obj.location, edge_dir );
        ap = obj.location + lambda * edge_dir;

    end

else
    
    %   Special case where source & target coincide.
    
    % Project position on edge.
    lambda = dot( source_pos - obj.location, edge_dir );
    ap = obj.location + lambda * edge_dir;
        
end

assert( any( ~isnan( ap ) ) )

end