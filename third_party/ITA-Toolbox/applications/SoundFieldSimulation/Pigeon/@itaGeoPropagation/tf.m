function [ freq_data_linear, valid ] = tf( obj, pp )
%TFS Calculates the transfer functions (tfs) of the (geometrical) propagation paths in frequency domain

if ~isfield( pp, 'propagation_anchors' )
    error( 'The propagation_path argument does not contain a field "propagation_anchors"' )
end

N = numel( pp.propagation_anchors );
if N < 2
    error( 'Propagation path has less than two anchor points, cannot calculate a transfer function' )
end

paths_distance = ita_propagation_path_length( pp );
if paths_distance / obj.c >  2 * obj.num_bins / obj.fs
    error( 'Propagation path length too long (%.1fm vs %.d bins), increase number of bins to generate transfer function for this propagation path', paths_distance, obj.num_bins )
end

freq_data_linear = obj.tf_atmospheric_absorption( paths_distance );

% Reflection & diffraction order
[ ro, do ] = ita_propagation_path_orders( pp );

valid =  true;
if obj.sim_prop.orders.combined > 0 && ( ro + do ) > obj.sim_prop.orders.combined
    %warning( 'Setting path %s as invalid, because reflection order %i plus diffraction order %i is exceeding configured combined order of %i', pp.identifier, ro, do, obj.sim_prop.orders.combined )
    valid = false;
    return
end

incident_spreading_loss_applied = false;

for n = 1 : N
    
    if isa( pp.propagation_anchors, 'cell' )
        anchor = pp.propagation_anchors{ n };
    else
        anchor = pp.propagation_anchors( n );
    end
    assert( strcmpi( anchor.class, 'propagation_anchor' ) )
    
    assert( isfield( anchor, 'anchor_type' ) )
    switch( anchor.anchor_type )
        
        case { 'source', 'emitter', 'receiver', 'sensor' }
            
            if n == N
                if isa( pp.propagation_anchors, 'cell' )
                    target_pos = pp.propagation_anchors{ n - 1 }.interaction_point;
                else
                    target_pos = pp.propagation_anchors( n - 1 ).interaction_point;
                end
                
                % If not already applied by a corresponding anchor type
                % (i.e. diffraction), include incident field now
                if ~incident_spreading_loss_applied
                    if ~obj.sim_prop.diffraction_enabled
                        effective_source_distance = distance_p; % whole distance in this case
                    else
                        effective_source_distance = ita_propagation_effective_source_distance( pp, n );
                    
                    end
                    phase_by_delay = obj.phase_delay( effective_source_distance );
                    spreading_loss = ita_propagation_spreading_loss( effective_source_distance, 'spherical' );
                    freq_data_linear = freq_data_linear .* phase_by_delay .* spreading_loss;
                    incident_spreading_loss_applied = true;
                end
                
            else
                
                if isa( pp.propagation_anchors, 'cell' )
                    target_pos = pp.propagation_anchors{ n + 1 }.interaction_point;
                else
                    target_pos = pp.propagation_anchors( n + 1 ).interaction_point;
                end
                
                % Check if sound power is set
                if n == 1
                    
                    p_factor = 1;
                    if isfield( anchor, 'sound_power' )
                        
                        r = 1;
                        rho_0 = 1.292; % Density of air
                        Z_0 = ( rho_0 * obj.c );
                        A = ( 4 * pi * r^2 );
                        I = anchor.sound_power / A;
                        p_factor = sqrt( I * Z_0 ); % Pressure factor @ 1m reference distance
                        
                    end
                    
                    assert( numel( p_factor ) == 1 ) % signal value scalar
                    freq_data_linear = freq_data_linear .* p_factor; % Apply factor corresponding to given sound power
                    
                end
                
            end
            
            target_position_relative = target_pos( 1:3 ) - anchor.interaction_point( 1:3 ); % Icoming or outgoing direction vector
            
            if obj.sim_prop.directivity_enabled
                freq_data_linear = freq_data_linear .* obj.tf_directivity( anchor, target_position_relative / norm( target_position_relative ) );
            end
            
        case 'specular_reflection'
            
            if n == 1 || n == N
                error( 'Detected a specular reflection at beginning or end of propagation path.' )
            end
            
            source_pos = pp.propagation_anchors{ n - 1 }.interaction_point;
            target_pos = pp.propagation_anchors{ n + 1 }.interaction_point;
            
            effective_source_position =  anchor.interaction_point - source_pos;
            target_position_relative =  target_pos - anchor.interaction_point;
            
            incident_direction_vec = effective_source_position / norm( effective_source_position );
            emitting_direction_vec = target_position_relative / norm( target_position_relative );
            
            if obj.sim_prop.reflection_enabled
                freq_data_linear = freq_data_linear .* obj.tf_reflection( anchor, incident_direction_vec, emitting_direction_vec );
            end
            
        case { 'outer_edge_diffraction', 'inner_edge_diffraction' }
            
            if n == 1 || n == N
                error( 'Detected a diffraction at beginning or end of propagation path.' )
            end
            
            % Diffraction values are summed up using the source /
            % receiver effective distances.
            % Spreading loss is spherical until first diffraction, after
            % that its a combination of spherical and cylindrical
            
            source_pos = pp.propagation_anchors{ n - 1 }.interaction_point( 1:3 );
            target_pos = pp.propagation_anchors{ n + 1 }.interaction_point( 1:3 );
                        
            source_direction = ( source_pos - anchor.interaction_point( 1:3 ) ) / norm( source_pos - anchor.interaction_point( 1:3 ) );
            target_direction = ( target_pos - anchor.interaction_point( 1:3 ) ) / norm( target_pos - anchor.interaction_point( 1:3 ) );
            
            if any( isnan( source_direction ) ) || any( isnan( target_direction ) )
                warning( 'Invalid path detected: %s (NaNs in edge diffraction calculation)', pp.identifier );
                valid = false;
                break;
            end
            
            effective_source_distance = ita_propagation_effective_source_distance( pp, n );
            effective_target_distance = ita_propagation_effective_target_distance( pp, n );
            effective_source_position = anchor.interaction_point( 1:3 ) + source_direction * effective_source_distance;
            effective_target_position = anchor.interaction_point( 1:3 ) + target_direction * effective_target_distance;
            
            assert( effective_source_distance > 0 && effective_target_distance > 0 )
            
            if obj.sim_prop.diffraction_enabled

                if ~obj.validate_diffraction( anchor, effective_source_position, effective_target_position )
                    warning( 'Invalid diffraction situation detected: %s (source or receiver inside wedge)', pp.identifier );
                    valid = false;
                    break;
                end

                if ~incident_spreading_loss_applied

                    % Propagation parameters
                    phase_by_delay = obj.phase_delay( effective_source_distance + effective_target_distance ); % Phase along full distance
                    spreading_loss = ita_propagation_spreading_loss( effective_source_distance, 'spherical' );
                    
                    % Apply fade-out on last diffraction
                    fadout_last_diffr = 1;
                    if do == obj.sim_prop.orders.combined
                        fadout_last_diffr = obj.fading_diffraction_utd( anchor, effective_source_position, effective_target_position ); % flip sign to cover phase inversion
                    end
                    
                    freq_data_linear = freq_data_linear .* phase_by_delay .* spreading_loss .* fadout_last_diffr;
                    
                    incident_spreading_loss_applied = true;
                    
                end
                
                [ diffr, valid ] = obj.tf_diffraction( anchor, effective_source_position, effective_target_position, obj.diffraction_model );
                                       
                if valid
                    freq_data_linear = freq_data_linear .* diffr;
                else
                    warning( 'Invalid diffraction calculation: %s', pp.identifier );
                    break
                end
                
            end
            
        otherwise
            
            sprintf( 'Detected unrecognized anchor type "%s", attempting to continue', anchor.anchor_type )
            
    end
    
end

if valid
    assert( incident_spreading_loss_applied )
else
    fprintf( 'Invalid path detected: %s\n', pp.identifier );
end

end
