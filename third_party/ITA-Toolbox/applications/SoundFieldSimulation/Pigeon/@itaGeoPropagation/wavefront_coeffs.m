function [ spreading_loss, diffraction_gain, delay, air_attenuation, object_interaction, source_wf_normal, receiver_wf_normal, valid, directivity ] = wavefront_coeffs( obj, pps )
%WAVEFRONT_COEFFS Derives acoustic parameters for given propagation paths.
%
%   Input:
%   pps                 Array of propagation path structs
%
%   Outputs (N=number of paths, M=number of bins):
%   spreading_loss      Spreading loss factors for each path (Nx1)
%   diffraction_gain    Gain compensating missing diffraction filter phase (Nx1)
%   delay               Propagation delay for each path (Nx1)
%   air_attenuation     Air attenuation frequency magnitudes for each path (MxN)
%   object_interaction  Object interaction (reflection + diffraction) frequency magnitudes (MxN)
%   source_wf_normal    Wavefront normal at source (Nx3)
%   receiver_wf_normal  Wavefront normal at receiver (Nx3)
%   valid               Boolean vector showing which propagation path is valid (Nx1)
%   directivity         [NOT_IMPLEMENTED] Source directivity frequency magnitudes (MxN)

%% Input checks
if ~isfield( pps, 'propagation_anchors' )
    error( 'The propagation_path argument does not contain a field "propagation_anchors"' )
end
assert( iscolumn(obj.freq_vector), 'Provide a valid frequency vector (must be a column vector)' )

%% Init
nPaths = numel(pps);
nFreqs = numel(obj.freq_vector);

spreading_loss = zeros(nPaths, 1);
diffraction_gain = zeros(nPaths, 1);
delay = zeros(nPaths, 1);
air_attenuation = zeros(nFreqs, nPaths);
object_interaction = zeros(nFreqs, nPaths);
source_wf_normal = zeros(nPaths, 3);
receiver_wf_normal = zeros(nPaths, 3);
valid = zeros(nPaths, 1);
directivity = []; %TODO

%% Calculation of parameters
for idxPath = 1:nPaths
    [ spreading_loss(idxPath), diffraction_gain(idxPath), delay(idxPath), air_attenuation(:, idxPath), object_interaction(:, idxPath), ...
        source_wf_normal(idxPath,:), receiver_wf_normal(idxPath,:), valid(idxPath) ] = ...
        wavefront_coeffs_single_path( obj, pps(idxPath) );
end


%% Calculation for single path
function [ spreading_loss, diffraction_gain, delay, air_attenuation, object_interaction, source_wf_normal, receiver_wf_normal, valid ] = wavefront_coeffs_single_path( obj, pp )

N = numel( pp.propagation_anchors );
if N < 2
    error( 'Propagation path has less than two anchor points, cannot calculate a transfer function' )
end

paths_distance = ita_propagation_path_length( pp );
delay = paths_distance  / obj.c;

freqs = obj.freq_vector;
spreading_loss = 1;
diffraction_gain = 1;
source_wf_normal = [ 0 0 0 ];
receiver_wf_normal = [ 0 0 0 ];

air_attenuation = ( 1 - ita_atmospheric_absorption_factor( freqs, paths_distance ) ); % Air absorption
object_interaction = ones(size(air_attenuation));

% Reflection & diffraction order
[ ro, do ] = ita_propagation_path_orders( pp );

valid = true;
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
                
                receiver_wf_normal = anchor.interaction_point(1:3) - target_pos(1:3);
                receiver_wf_normal = receiver_wf_normal' / norm(receiver_wf_normal);
                
                % If not already applied by a corresponding anchor type
                % (i.e. diffraction), include incident field now
                if ~incident_spreading_loss_applied
                    if ~obj.sim_prop.diffraction_enabled
                        effective_source_distance = ita_propagation_path_length( pp ); % whole distance in this case
                    else
                        effective_source_distance = ita_propagation_effective_source_distance( pp, n );
                    end
                    spreading_loss = spreading_loss * ita_propagation_spreading_loss( effective_source_distance, 'spherical' );
                    incident_spreading_loss_applied = true;
                end
                
            else
                
                if isa( pp.propagation_anchors, 'cell' )
                    target_pos = pp.propagation_anchors{ n + 1 }.interaction_point;
                else
                    target_pos = pp.propagation_anchors( n + 1 ).interaction_point;
                end
                
                if n == 1                    
                    source_wf_normal = target_pos(1:3) - anchor.interaction_point(1:3);
                    source_wf_normal = source_wf_normal' / norm(source_wf_normal);
                end
                
            end
            
            target_position_relative = target_pos( 1:3 ) - anchor.interaction_point( 1:3 ); % Incoming or outgoing direction vector
            
            if obj.sim_prop.directivity_enabled && false
                directivity = obj.mags_directivity( anchor, target_position_relative / norm( target_position_relative ), freqs ); %  todo
            end
            
        case 'specular_reflection'
            
            if n == 1 || n == N
                error( 'Detected a specular reflection at beginning or end of propagation path.' )
            end
            
            source_pos = pp.propagation_anchors{ n - 1 }.interaction_point;
            target_pos = pp.propagation_anchors{ n + 1 }.interaction_point;
            
            effective_source_position =  anchor.interaction_point(1:3) - source_pos(1:3);
            target_position_relative =  target_pos(1:3) - anchor.interaction_point(1:3);
            
            incident_direction_vec = effective_source_position / norm( effective_source_position );
            emitting_direction_vec = target_position_relative / norm( target_position_relative );
            
            if obj.sim_prop.reflection_enabled && false
                object_interaction = object_interaction .* obj.mags_reflection( anchor, incident_direction_vec, emitting_direction_vec, freqs ); % todo
                %object_interaction = object_interaction .* db2mag( -0.9 );
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
            
            if obj.sim_prop.diffraction_enabled
                
                if ~obj.validate_diffraction( anchor, effective_source_position, effective_target_position )
                    warning( 'Invalid diffraction situation detected: %s (source or receiver inside wedge)', pp.identifier );
                    valid = false;
                    break;
                end
                
                if ~incident_spreading_loss_applied
                    
                    spreading_loss = spreading_loss *  ita_propagation_spreading_loss( effective_source_distance, 'spherical' );
                    diffr_sign = obj.sign_diffraction_utd( anchor, effective_source_position, effective_target_position ); % flip sign to cover phase inversion
                    diffraction_gain = diffraction_gain * diffr_sign;
                    
                    % Apply fade-out on last diffraction
                    if do == obj.sim_prop.orders.combined
                        diffr_fading = obj.fading_diffraction_utd( anchor, effective_source_position, effective_target_position ); % flip sign to cover phase inversion
                        %ita_diffraction_utd_illumination_fadeout( w, s, r3, false ); % fade out?
                        diffraction_gain = diffraction_gain * diffr_fading;
                    end

                    incident_spreading_loss_applied = true;
                    
                end
                
                [ temp_diffraction_mags, valid ] = obj.mags_diffraction_utd( anchor, effective_source_position, effective_target_position, freqs );
                                                
                if valid
                    if any( isnan( temp_diffraction_mags ) )
                        error( 'Detected NaNs in wavefront coefficint calculation during edge diffraction, segmend %i of path %s', n, pp.identifier );
                    end
                    object_interaction = object_interaction .* temp_diffraction_mags;
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
    if any( isnan( object_interaction ) )
        error( 'Detected NaNs in wavefront coefficint calculation of path %s', pp.identifier );
    end
else
    fprintf( 'Invalid path detected: %s\n', pp.identifier );
end
