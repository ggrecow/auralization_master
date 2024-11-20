%% PIGEON URBAN SIMULATION: AURALIZATION EXAMPLE
% Example how to use the pigeon interface itaPigeonProject for
% the auralization of a dynamic scene. It utilizies the
% OutdoorNoise renderer of the Virtual Acoustics auralization frame
% work
% 
% The auraliazion requires:
% 2. pigeon interface (>= v2021a) for the urban simulation
% 3. Virtual Acoustics (>= v2021a) inlcuding BinauralOutdoorNoise renderer

%% Settings
%CAU simulation settings
forceCalculation = false;
simResultFolder = fullfile(cd,'pigeon_results');
freqVector = ita_ANSI_center_frequencies;
urbanGeometryFile = which('pigeon_test.skp');

%VA / audio settings
rendererName = 'MyBinauralOutdoorUrban';
conf_path = which( 'VACore.OutdoorNoise.pigeon.recording.ini' );
output_filename = 'pigeon_auralization_example.wav';
core_config = IniConfig();
core_config.ReadFile(conf_path);
samplingRate = core_config.GetValues('Audio driver','Samplerate');
blockSize = core_config.GetValues('Audio driver','Buffersize');

%Receiver position
if true
    receiverPosition = [7 0.8 1.7];
    receiverView = [ -1 0 0];
    receiverUp = [0 0 1];
else
    receiverPosition = [1 7 1.7];
    receiverView = [ 1 0 0];
    receiverUp = [0 0 1];
end

%Car trajectory
carVelocity = 30 / 3.6; %30 km/h -> m/s
carStartPos = [3 16 1.7];
carEndPos = [3 0 1.7];
tMax = norm(carEndPos - carStartPos) / carVelocity; %[s]

%Sound source power
% soundPowerLevel = 95; %[dB]
% soundPower = 1e-12 * 10^(soundPowerLevel/10); %[W]
soundPower = 20e-3;

%% Dependent settings
%VA / audio
deltaT = blockSize/samplingRate;
updateRate = 1 / deltaT;

%Car trajectory
t = (0:deltaT:tMax)';
carPositions = [...
    interp1([0,tMax], [carStartPos(1), carEndPos(1)], t),...
    interp1([0,tMax], [carStartPos(2), carEndPos(2)], t),...
    interp1([0,tMax], [carStartPos(3), carEndPos(3)], t)...
    ];

%% Run simulation if not done and write into file

pigeon = itaPigeonProject;
pigeon.geometry_file_path = urbanGeometryFile;
propagation = itaGeoPropagation;

sim_result_files = pigeon.run_dynamic_scenario(carPositions, receiverPosition, simResultFolder, forceCalculation);
auralization_frame_files = propagation.calculate_scene_auralization(pigeon, sim_result_files, simResultFolder, freqVector, forceCalculation);

nFiles = numel(auralization_frame_files);

%% Create VA object (opens GUI to select VA path if required)
va = VA();

%% Start VA server
[ basename, basepath ] = uigetfile('VAServer.exe','Select VAServer executable','VAServer');

if basename == 0
    error('No VAServer executable selected')
end

va_basepath = fullfile( basepath, '..' );

va_args = [ 'localhost:12340 "' conf_path '"'];
os_call = [ basepath basename ' ' va_args ' rc &' ];
% ... starts in remote control mode, use va.shutdown_server to stop server and export WAV file
system(os_call);

%% Initialize VA scene
va.connect( 'localhost' );

va.set_output_gain( 0.5 )

tmp_filename = 'pigeon_auralization_tmp.wav';
va.set_rendering_module_parameters( rendererName, struct('RecordOutputFileName', tmp_filename) );

va.add_search_path( fullfile(va_basepath,'data') );
va.add_search_path( cd );

%Signal source
sourceSignal = va.create_signal_source_buffer_from_file( 'CarNoise.wav' );
va.set_signal_source_buffer_looping( sourceSignal, true );
va.set_signal_source_buffer_playback_action( sourceSignal, 'play' )

%Sound source
S = va.create_sound_source( 'itaVA_Source' );
va.set_sound_source_signal_source( S, sourceSignal );
va.set_sound_source_sound_power( S, soundPower);
va.set_sound_source_orientation_view_up(S,[0 0 -1],[0 1 0]);

%Receiver
R = va.create_sound_receiver( 'itaVA_Listener' );
va.set_sound_receiver_position( R, ita_matlab2openGL( receiverPosition ) )
va.set_sound_receiver_orientation_view_up( R, ita_matlab2openGL(receiverView), ita_matlab2openGL(receiverUp));
H = va.create_directivity_from_file( '$(DefaultHRIR)' );
va.set_sound_receiver_directivity( R, H );

%% Preprocessing so that VDL can fill up
manual_clock = 0;
va.set_core_clock( 0 );
tPreProc = 2;
nPreBlocks = ceil( tPreProc * samplingRate / blockSize );
nPreSamples = nPreBlocks*blockSize;

for idFrame = 1:nPreBlocks
    va.call_module( 'virtualaudiodevice', struct( 'trigger', true ) );
    
    manual_clock = manual_clock + deltaT;
    va.call_module( 'manualclock', struct( 'time', manual_clock ) );
end

%% Run actual auralization

wb = itaWaitbar(nFiles,'Running Auralization...');
for idFrame = 1:nFiles
   
    % Make source_pos, receiver_pos and paths_update available for this
    % frame
    load( auralization_frame_files{idFrame} , 'paths_update', '-mat' );
    path_names = fieldnames(paths_update);
          
    % Update all propagation paths
    source_pos_OpenGL = ita_matlab2openGL( carPositions(idFrame, :) );
    receiver_pos_OpenGL = ita_matlab2openGL(receiverPosition);
    
    va.set_sound_source_position( S, source_pos_OpenGL );
    va.set_sound_receiver_position( R, receiver_pos_OpenGL );
    
    va.set_rendering_module_parameters( rendererName, paths_update );
    
    % Process audio chain by incrementing one block
    va.call_module( 'virtualaudiodevice', struct( 'trigger', true ) );
       
    manual_clock = manual_clock + deltaT;
    va.call_module( 'manualclock', struct( 'time', manual_clock ) );

    if isvalid(wb); wb.inc(); end
end
if isvalid(wb); wb.close(); end

pause(1)
va.shutdown_server();
pause(5); %Wait for VA to write file
va.disconnect();

%% Remove time used to fill up VDL
audioSamples = audioread(tmp_filename);
audioSamples = audioSamples(nPreSamples+1:end, :);
audiowrite(output_filename, audioSamples, samplingRate, 'BitsPerSample', 32);
delete(tmp_filename);