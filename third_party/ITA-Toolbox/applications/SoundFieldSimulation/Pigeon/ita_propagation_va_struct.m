function paths_update = ita_propagation_va_struct(pps, pps_old, S, R, spreading_loss, diffraction_gain, delay, source_wf_normals, receiver_wf_normals, air_attenuation, object_interaction, directivity, b2OpenGL)
%ITA_PROPAGATION_VA_STRUCT Creates a struct based on the acoustic
%parameters derived from a set of propagation paths, which can be sent to
%the OutdoorNoiseRenderer or BinauralOutdoorNoiseRenderer of Virtual Acoustics (VA).
%
%   Inputs (default):
%   pps                 Struct with propagation paths
%   pps_old             Struct with propagation paths of last time frame (leave empty on first call)
%   S                   VA source ID [int]
%   R                   VA receiver ID [int]
%   spreading_loss      Gain factor representing spreading loss
%   diffraction_gain    Gain factor compensating missing phase in diffraction filter
%   delay               Propagation delay
%   source_wf_normals   Wavefront normals at source [Nx3 vector]
%   receiver_wf_normals Wavefront normals at receiver [Nx3 vector]
%   air_attenuation     Third octave frequency magnitudes for air attenuation [31xN vector]
%   object_interaction  Third octave frequency magnitudes for object interaction (reflection + diffraction) [31xN vector]
%   directivity ([])    Third octave frequency magnitudes for source directivity [31xN vector]
%   b2OpenGL (false)    If true, converts all geometric vectors to OpenGL coordinates

if nargin < 12
    directivity = [];
end
if nargin < 13
    b2OpenGL = true;
end
nPaths = numel(pps);
paths_update = struct();

%% Input checks
assert( isstruct(pps) && isfield(pps, 'propagation_anchors'), 'pps must be a propagation paths struct!' )
assert( isempty( pps_old ) || ( isstruct(pps_old) && isfield(pps_old, 'propagation_anchors') ), 'pps_old must be a propagation paths struct!' )
assert( isnumeric(S) && isscalar(S) && mod(S,1) == 0, 'S must be an integer!' )
assert( isnumeric(R) && isscalar(R) && mod(R,1) == 0, 'R must be an integer!' )
assert( isnumeric(air_attenuation) && size(air_attenuation, 1) == 31, 'Air attenuation must be a 31xN matrix.' )
assert( isnumeric(object_interaction) && size(object_interaction, 1) == 31, 'Object interaction must be a 31xN matrix.' )
assert( isnumeric(source_wf_normals) && size(source_wf_normals, 2) == 3, 'Wavefront normals must be a Nx3 matrix.' )
assert( isnumeric(receiver_wf_normals) && size(receiver_wf_normals, 2) == 3, 'Wavefront normals must be a Nx3 matrix.' )
assert( isnumeric(spreading_loss) && numel(spreading_loss) == nPaths, 'Number of spreading loss values must match number of paths.' )
assert( isnumeric(diffraction_gain) && numel(diffraction_gain) == nPaths, 'Number of gain values must match number of paths.' )
assert( isnumeric(delay) && numel(delay) == nPaths, 'Number of delay values must match number of paths.' )
assert( size(air_attenuation, 2) == nPaths, 'Number of air absorption magnitude vectors must match number of paths.' )
assert( size(object_interaction, 2) == nPaths, 'Number of object interaction magnitude vectors must match number of paths.' )
assert( size(source_wf_normals, 1) == nPaths, 'Number of wavefront normals must match number of paths.' )
assert( size(receiver_wf_normals, 1) == nPaths, 'Number of wavefront normals must match number of paths.' )
assert( isempty(directivity) || isnumeric(directivity) && size(directivity) == size(air_attenuation), 'Directivity must be empty or a 31xN matrix.' )

%% Add path identifiers if required
if ~isfield( pps, 'identifier' )
    pps = ita_propagation_paths_add_identifiers( pps );
end
if ~isempty( pps_old ) && ~isfield( pps_old, 'identifier' )
    pps_old = ita_propagation_paths_add_identifiers( pps_old );
end

%% Get deleted paths
if isempty(pps_old)
    pps_del = [];
else
    [ ~, pps_del, ~ ] = ita_propagation_paths_diff( pps_old, pps );
end

%% Delete non-available paths
for iPath = 1:numel( pps_del )
    pu = struct(); % Path update
    pu.source = S;
    pu.receiver = R;
    pu.identifier = pps_del( iPath ).identifier;
    pu.delete = true;
    paths_update.( strcat( 'path_', pu.identifier ) ) = pu;  
end

%% Valid paths
path_directivity = [];
for idx = 1:numel(pps)
    if ~isempty(directivity); path_directivity = directivity(:, idx); end
    pu = single_path_va_struct( pps(idx), S, R, spreading_loss(idx), diffraction_gain(idx),...
          delay(idx), source_wf_normals(idx,:), receiver_wf_normals(idx,:),...
         air_attenuation(:, idx), object_interaction(:, idx), path_directivity, b2OpenGL );
    paths_update.( strcat( 'path_', pu.identifier ) ) = pu;
end

%% Add source and receiver ID
paths_update.sound_source_id = S;
paths_update.sound_receiver_id = R;


function pu = single_path_va_struct(pp, S, R, spreading_loss, diffraction_gain, delay, source_wf_normal, receiver_wf_normal, air_attenuation, object_interaction, directivity, b2OpenGL)

%% Input checks
assert( ~isempty(pp.identifier), 'All path structs must have a valid identifier!' )


%% Extract last interaction point
if isa(pp.propagation_anchors, 'struct')
    lastIntPoint = pp.propagation_anchors( end-1 ).interaction_point( 1:3 )'; % next to last anchor
else
    lastIntPoint = pp.propagation_anchors{ end-1 }.interaction_point( 1:3 )'; % next to last anchor
end

%% Coordinate conversion
if b2OpenGL
    lastIntPoint = ita_matlab2openGL( lastIntPoint );
    source_wf_normal = ita_matlab2openGL( source_wf_normal );
    receiver_wf_normal = ita_matlab2openGL( receiver_wf_normal );
end

%% Create VA struct
pu.sound_source_id = S;
pu.sound_receiver_id = R;
pu.identifier = pp.identifier;

pu.spreading_loss = spreading_loss;
pu.diffraction_gain = diffraction_gain;
pu.propagation_delay = delay;

pu.air_attenuation_third_octaves = air_attenuation;
pu.object_interaction_third_octaves = object_interaction;
if(~isempty(directivity))
    pu.directivity_third_octaves = directivity;
end

[ reflection_order, diffraction_order ] = ita_propagation_path_orders( pp );
pu.reflection_order = reflection_order;
pu.diffraction_order = diffraction_order;

pu.source_wavefront_normal = source_wf_normal;
pu.receiver_wavefront_normal = receiver_wf_normal;

pu.delete = false;
pu.audible = spreading_loss ~= 0;

%% Backwards compatibility (BinauralOutdoorNoiseRenderer)

total_mags = air_attenuation .* object_interaction;
if(~isempty(directivity))
    total_mags = total_mags .* directivity;
end

pu.source = S;
pu.receiver = R;

pu.gain = spreading_loss * diffraction_gain;
pu.delay = delay;
pu.frequency_magnitudes = total_mags;

pu.position = lastIntPoint;