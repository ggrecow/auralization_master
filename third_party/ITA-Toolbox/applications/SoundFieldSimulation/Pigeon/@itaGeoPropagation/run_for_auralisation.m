function [ paths_update ] = run_for_auralisation( obj, b2OpenGL )
%RUN_FOR_AURALIZTAION Calculates the auralisation coefficients of each
%(geometrical) propagation path and creates a struct that can be sent to
%Virtual Acoustics (VA).
%
%   Inputs (default):
%   b2OpenGL (false)    If true, converts all geometric vectors to OpenGL coordinates

%% Checks
assert( isfield(obj.pps, 'propagation_anchors'), 'Propagation paths struct (pps) must be set first!' )
assert( isempty(obj.pps_old) || isfield(obj.pps_old, 'propagation_anchors'), 'Property "pps_old" must either be empty or a propagation paths struct with paths of last iteration.' )

%% Default inputs
if nargin < 2; b2OpenGL = false; end

if ~isfield( obj.pps, 'identifier' )
    obj.pps = ita_propagation_paths_add_identifiers( obj.pps );
end
if ~isfield( obj.pps_old, 'identifier' )
    obj.pps_old = ita_propagation_paths_add_identifiers( obj.pps_old );
end

%% Limit combined reflection / diffraction order
if obj.sim_prop.orders.combined >= 0
    pps = ita_propagation_path_filter_combined_order( obj.pps, obj.sim_prop.orders.combined );
    pps_old = ita_propagation_path_filter_combined_order( obj.pps_old, obj.sim_prop.orders.combined );
else
    pps = obj.pps;
    pps_old = obj.pps_old;
end

%% Process
[ spreading_loss, diffraction_gain, delay, air_attenuation, object_interaction, source_wf_normals, receiver_wf_normals, valid ] = obj.wavefront_coeffs( pps );
%TODO: Use the valid boolean vector to filter out invalid sound paths
paths_update = ita_propagation_va_struct(pps, pps_old, obj.source_id, obj.receiver_id, spreading_loss, diffraction_gain, delay, source_wf_normals, receiver_wf_normals, air_attenuation, object_interaction, [], b2OpenGL);
