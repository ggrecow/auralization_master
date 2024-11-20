%% COMBINED ATMOSPHERIC URBAN SIMULATION: AURALIZATION EXAMPLE
% Example how to use the combined simulation interface itaCauSimulation for
% the auralization of a dynamic scene. It utilizies the
% BinauralOutdoorNoise renderer of the Virtual Acoustics auralization frame
% work
% 
% The auraliazion requires:
% 1. ARTMatlab (>= v2021b) framework for the atmospheric simulation
% 2. pigeon interface (>= v2021a) for the urban simulation
% 3. Virtual Acoustics (>= v2021a) inlcuding BinauralOutdoorNoise renderer

%% Settings
%CAU simulation settings
forceCalculation = false;
simResultFolder = fullfile(cd,'cau_results');

freqVector = ita_ANSI_center_frequencies;
urbanGeometryFile = which('streetCanyon.skp');

%VA / audio settings
samplingRate = 44100;
blockSize = 1024;
rendererName = 'MyBinauralOutdoorNoise';
conf_path = which( 'VACore.OutdoorNoise.CAU.recording.ini' );
% conf_path = which( 'VACore.OutdoorNoise.CAU.recording.v2022a.ini' );
output_folder = cd;
output_filename = fullfile(output_folder, 'cau_auralization_example.wav');

%Receiver position
receiverPosition = [17.5,100,1.8];
receiverView = [0 1 0];
receiverUp = [0 0 1];

%Aircraft trajectory
aircraftVelocity = 350 / 3.6; %350 km/h -> m/s
aircraftAltitude = 1000; %[m]
aircraftFlightDir = [1 0 0];
trajectoryCenterPos = [receiverPosition(1:2) aircraftAltitude];
tMax = 15; %[s]

%Sound source power
soundPower = 20; %[W]

%% Dependent settings
%VA / audio
deltaT = blockSize/samplingRate;
updateRate = 1 / deltaT;

%Aircraft trajectory
t = (0:deltaT:tMax)';
aircraftMovement = aircraftVelocity * t * aircraftFlightDir;
aircraftStartPos = trajectoryCenterPos - aircraftMovement(end, :)/2;
aircraftPositions = aircraftStartPos + aircraftMovement;

%% Run simulation if not done and write into file

cauSimulation = itaCauSimulation;
cauSimulation.urbanPropagation.MaxDiffractionOrder = 2;
cauSimulation.urbanPropagation.MaxReflectionOrder = 2;
cauSimulation.urbanPropagation.MaxCombinedOrder = 2;
cauSimulation.urbanPropagation.geometry_file_path = urbanGeometryFile;
cauProp = itaCauPropagation( cauSimulation );

sim_result_files = cauSimulation.DynamicScenario(aircraftPositions, receiverPosition, updateRate, simResultFolder, forceCalculation);
auralization_frame_files = cauProp.AuralizationScenarioParameters(sim_result_files, simResultFolder, freqVector, forceCalculation);

nFiles = numel(auralization_frame_files);

%% Create VA object (opens GUI to select VA path if required)
va = VA();

%% Start VA server
[ basepath, basename, ext ]= fileparts( which( 'VAServer.exe') );
[ va_basepath, ~, ~ ]= fileparts( basepath );

va_args = [ 'localhost:12340 "' conf_path '"'];
os_call = [ which( 'VAServer.exe' ) ' ' va_args ' rc &' ];
% ... starts in remote control mode, use va.shutdown_server to stop server and export WAV file
system(os_call);

%% Initialize VA scene
va.connect( 'localhost' );

tmp_filename = 'cau_auralization_tmp.wav';
va.set_rendering_module_parameters( rendererName, struct('RecordOutputBaseFolder', output_folder, 'RecordOutputFileName', tmp_filename) );

va.add_search_path( fullfile(va_basepath,'data') );

%Signal source
sspt_jet_engine_conf.class = 'jet_engine';
sourceSignal = va.create_signal_source_prototype_from_parameters( sspt_jet_engine_conf );
va.set_signal_source_parameters( sourceSignal, struct( 'rpm', 1300 ) )

%Sound source
S = va.create_sound_source( 'itaVA_Source' );
va.set_sound_source_signal_source( S, sourceSignal );
va.set_sound_source_sound_power( S, soundPower);
va.set_sound_source_orientation_view_up(S,[0 0 -1],[0 1 0]);

%Receiver
R = va.create_sound_receiver( 'itaVA_Listener' );
va.set_sound_receiver_position( R, receiverPosition )
va.set_sound_receiver_orientation_view_up( R, receiverView, receiverUp);
H = va.create_directivity_from_file( '$(DefaultHRIR)' );
va.set_sound_receiver_directivity( R, H );

%% Preprocessing so that VDL can fill up
manual_clock = 0;
va.set_core_clock( 0 );
tPreProc = 10;
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
          
    % Update all propagation paths
    source_pos_OpenGL = aircraftPositions(idFrame, :);
    receiver_pos_OpenGL = receiverPosition;
    
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