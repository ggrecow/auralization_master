function [ valid_pps, invalid_pps ] = ita_propagation_paths_filter_valid( pps, geo_eps )

if ~isfield( pps, 'propagation_anchors' ) % not a list but only one path
    error( 'Need a propagation path or path list' )
end

if nargin < 2
    geo_eps = 1e-6;
end

N = numel( pps );
valid_paths = 1:N;

for n = 1:N

    propagation_path = pps( n );
    
    M = numel( propagation_path.propagation_anchors );
    if M < 2
        error( 'Propagation path has less than two anchor points, cannot calculate a transfer function' )
    end
    
    
    for m = 1 : M-1
        
        if isa( propagation_path.propagation_anchors, 'cell' )
            cur_anchor_pos = propagation_path.propagation_anchors{ m }.interaction_point( 1:3 );
            next_anchor_pos = propagation_path.propagation_anchors{ m + 1 }.interaction_point( 1:3 );
            if m > 1
                prev_anchor_pos = propagation_path.propagation_anchors{ m - 1 }.interaction_point( 1:3 );
            end
            cur_anchor = propagation_path.propagation_anchors{ m };
        else
            cur_anchor_pos = propagation_path.propagation_anchors( m ).interaction_point( 1:3 );
            next_anchor_pos = propagation_path.propagation_anchors( m + 1 ).interaction_point( 1:3 );
            if m > 1
                prev_anchor_pos = propagation_path.propagation_anchors( m - 1 ).interaction_point( 1:3 );
            end
            cur_anchor = propagation_path.propagation_anchors( m );
        end
        
        assert( isfield( cur_anchor, 'anchor_type' ) )
        if m == 1
            if ~( strcmpi( cur_anchor.anchor_type, 'emitter' ) || strcmpi( cur_anchor.anchor_type, 'source' ) )
                error( 'Removing path %i, first anchor was not of anchor_type ''emitter''', n )
            end
        elseif m == M
            if ~( strcmpi( cur_anchor.anchor_type, 'sensor' ) || strcmpi( cur_anchor.anchor_type, 'receiver' ) )
                error( 'Removing path %i, first anchor was not of anchor_type ''sensor''', n )
            end
        end
            
        % Far field condition
        far_field_distance =  norm( next_anchor_pos - cur_anchor_pos );
        if far_field_distance < geo_eps
            valid_paths( n ) = 0;
            warning( 'Removing path %i because anchor %i/%i was violating far field condition with next anchor point, distace was %.1dm', n, m, M, far_field_distance )
            break
        end
        
        % Diffraction source/target locations and apex point validation
        switch( cur_anchor.anchor_type )
        case { 'outer_edge_diffraction', 'inner_edge_diffraction' }
            
            if m == 1 || m == M
                error( 'Detected a diffraction at beginning or end of propagation path.' )
            end
            
            n1 = cur_anchor.main_wedge_face_normal( 1:3 );
            n2 = cur_anchor.opposite_wedge_face_normal( 1:3 );
            start_point = cur_anchor.vertex_start( 1:3 );
            end_point = cur_anchor.vertex_end( 1:3 );

            if abs( cross( n1, n2 ) ) < itaInfiniteWedge.set_get_geo_eps
                len = norm( end_point - start_point );
                edge_dir = ( end_point - start_point ) / len;
                w = itaSemiInfinitePlane( n1, start_point, edge_dir );
            else
                if strcmpi( cur_anchor.anchor_type, 'outer_edge_diffraction' )
                    w = itaInfiniteWedge( n1, n2, start_point );
                else
                    w = itaInfiniteWedge( n1, n2, end_point ); % flips edge direction
                end
            end
            
            if ~w.point_outside_wedge( prev_anchor_pos ) || ~w.point_outside_wedge( next_anchor_pos )
                valid_paths( n ) = 0;
                warning( 'Removing path %i with a diffraction wedge %i/%i, origin or target was located inside wedge', n, m, M )
                break
            end
            
            apex_point_approx = w.apex_point_approx( prev_anchor_pos, next_anchor_pos );
            d1 = norm( prev_anchor_pos - cur_anchor.interaction_point( 1:3 ) ) + norm( cur_anchor.interaction_point( 1:3 ) - next_anchor_pos );
            d2 = norm( prev_anchor_pos - apex_point_approx ) + norm( apex_point_approx - next_anchor_pos );                       
            distance_err1 = abs( d1 - d2 ) / abs( d1 + d2 );
            distance_cond1 = distance_err1 < 0.1;
            
            
            apex_point_far_field = w.get_aperture_point_far_field( prev_anchor_pos, next_anchor_pos );
            distance_err2 = norm( apex_point_far_field - cur_anchor.interaction_point( 1:3 ) );
            distance_cond2 = distance_err2 < 1.0;
            
            if ~( distance_cond1 || distance_cond2 )
                if ~distance_cond1
                    warning( 'Removing path %i because anchor %i/%i relative detour has a factor error of %.1fm', n, m, M, distance_err1 )
                end
                
                if ~distance_cond2
                    warning( 'Removing path %i because anchor %i/%i is violating far field conditions, error distance was %.1fm', n, m, M, distance_err2 )
                end
                
                valid_paths( n ) = 0;
                break                
            end
            
            % too similar source / receiver incidence (bouncing between parallel edges) are
            % toublesome in utd, remove.
            alpha_i = w.angle_main_face( prev_anchor_pos );
            alpha_d = w.angle_main_face( next_anchor_pos );
            if abs( alpha_d - alpha_i ) < deg2rad( 1 ) % degree
                warning( 'Removing path %i because detected a too steep in/out diffraction angle for source and receiver, UTD will produce high levels', n )
                valid_paths( n ) = 0;
                break
            end
            
            % Remove inner edges
            if strcmpi( cur_anchor.anchor_type, 'inner_edge_diffraction' )
                warning( 'Removing path %i because detected an inner edge diffraction anchor', n )
                valid_paths( n ) = 0;
                break
            end
                    
        case { 'source', 'emitter', 'receiver', 'sensor' }            
            
        case 'specular_reflection'
            
        otherwise
            
        end
    end
    
end

valid_paths_idx = valid_paths( valid_paths > 0 );
valid_pps = pps( valid_paths_idx );

if nargout > 1
    invalid_pps = pps( ~valid_paths );
end

end