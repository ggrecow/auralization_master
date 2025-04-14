function OutputAuralization = MASTER_auralization_EngineAirframe(input, flight_profile, smoothing, input_type, tag_auralization, show)
% function OutputAuralization = MASTER_auralization_EngineAirframe(input, flight_profile, smoothing, input_type, tag_auralization, show)
%
% INPUTS:
%
%   input : struct
%   struct containing all processed data per time step in   
%                   the format - data{nTimeSteps,nObserver} - contains
%                   db(Z) and dB(A) data. This is the same as the <data>
%                   output from the <PANAM_SQAT_data_conversion>. As a best practice, it is 
%                   recommended to use here the <data> struct which was
%                   trimmed by the <prepare_input_SQ> function 
%
%   flight_profile : struct
%   contain variables related to flight profile for propagation using ray
%   tracing (only used if input_type = 'emission')
%
%   smoothing : scalar
%   number of spectral smoothing iterations during the conversion from toc to
%   narrowband (can be ZERO)
%
%   input_type : string
%   controls the type of data used as input for auralization. Input values are:
%
%   - 'emission' - ray tracing will be used to compute sound propagation and
%     emission time will be used to deal with time vectors
%
%   - 'immission' - auralization computed directly from immission values
%      (i.e. propagation already included) and retarted time will be used to
%      compute the auralization time vector
%
%   tag_auralization : string
%   containg the main path (where files should be saved) and their corresponding basic info 
%   (i.e. emission/immission & aircraft & procedure & Receiver #) to be
%   included in the name of the saved files related to auralization processes
%
%   show : logical (boolean)
%   optional parameter for figures (results) display
%   'false' (disable, default value) or 'true' (enable)
%
% OUTPUTS
%
% Author: Gil Felix Greco, Braunschweig 20.06.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global input_file

% reference pressure [Pa]
global pref 
pref = 20e-6;         

% sampling frequency used for auralization
global fs

if isfield( input_file, 'sampling_freq' )
    fs = str2double( input_file.sampling_freq ); % get <fs> from <input_file>
else
    fs = 48000; %default value
end    

% define default for <binaural_signal>
if isfield( input_file, 'binaural_signal' )
    binaural_signal = str2double ( input_file.binaural_signal );   % boolean : 0= only direct path; 1 = direct path + 1st order reflection
else
    binaural_signal = 0; % default value
end

% time resolution from PANAM. Assumes it is constant (source noise prediction dt from PANAM is usually .5 seconds and cte) - ACHTUNG: this is the case only for sound source   
global dt_panam
dt_panam =  input{2,1}.source_time - input{1,1}.source_time;  

%% get time vector

[~, time_PANAM_auralization, time]  = getAuralizationTime( input, input_type );

%% Synthesis of tonal signal (engine - fan harmonics)

tag_source = 'fan_harmonics';

input_fan_harmonics = getTonalInput( input, time_PANAM_auralization, tag_source, tag_auralization );
TonalSignal_FanHarmonics = tonalSynthesis( input_fan_harmonics, time_PANAM_auralization, time, show, tag_auralization, tag_source, input_type ); 

clear input_fan_harmonics;

%% Synthesis of tonal signal (engine - buzzsaw)

tag_source='buzzsaw';

input_buzzsaw = getTonalInput( input, time_PANAM_auralization, tag_source, tag_auralization );
TonalSignal_Buzzsaw = tonalSynthesis( input_buzzsaw, time_PANAM_auralization, time, show, tag_auralization, tag_source, input_type ); 

clear input_buzzsaw;

%% Synthesis of broadband signal (engine)

tag_source = 'engine';
BroadbandSignal_engine = broadbandSynthesis_smooth( input, tag_source, time_PANAM_auralization, time, smoothing, show, tag_auralization, input_type );

%% engine signal ( BBN-engine + tonal-fan harmonics + tonal-buzzsaw)

auralizedEngineSignal = TonalSignal_FanHarmonics + TonalSignal_Buzzsaw + BroadbandSignal_engine;

%%%% check plot (spectrogram) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if show == true    
    tag_title =  [ 'OUTPUT - Spectrogram of auralized engine signal - ' input_type ' (BBN engine+fan_harmonics+buzzsaw)' ];
    tag_save =  '_engineSignal_Spectrogram' ;
    PLOT_spectrogram( auralizedEngineSignal, fs, tag_title, tag_auralization, tag_save) ;    
else
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Synthesis of broadband signal (airframe)

tag_source = 'airframe';
auralizedAirframeSignal = broadbandSynthesis_smooth( input, tag_source, time_PANAM_auralization, time, smoothing, show, tag_auralization, input_type );

%% Synthesis of broadband signal

auralizedOverallSignal = zeros( length( auralizedEngineSignal ),1 ); % pre allocate memory
auralizedOverallSignal = auralizedEngineSignal + auralizedAirframeSignal;

%%%% check plot (spectrogram) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if show == true    
    tag_title = [ 'OUTPUT - Spectrogram of auralized signal - ' input_type ' (engine + airframe)' ];
    tag_save =  '_overallSignal_spectrogram' ;
    PLOT_spectrogram(auralizedOverallSignal, fs, tag_title, tag_auralization, tag_save);    
else
end

% sound(EngineSignal,44100,24);
% sound(AirframeSignal,44100,24);
% sound(OverallSignal,44100,24);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Post process auralized signal
% - auralized emission signal:  apply propagation and save output (.wav and vector)

% case 'emission' % synthesized source signals are propagated, auralized and saved

% save final (EMISSION) auralized signal
OutputAuralization.emission.overallSignal = auralizedOverallSignal;

OutputAuralization.emission.engineSignal = auralizedEngineSignal;
OutputAuralization.emission.tonalSignalFanHarmonics = TonalSignal_FanHarmonics;
OutputAuralization.emission.tonalSignalBuzzsaw = TonalSignal_Buzzsaw;
OutputAuralization.emission.broadbandSignal_engine = BroadbandSignal_engine;

OutputAuralization.emission.auralizedAirframeSignal = auralizedAirframeSignal;

% plot propagation graphs?
show_propagation = 1; % plot propagation results

% FFT length used for freq vector of atmospheric transfer
% function and as nTaps of the FIR filter of atmospheric transfer function
% nfft = round(fs*dt_panam);
% nfft = str2double( input_file.nfft );
nfft = 2^14; % to yield df=3 Hz. needs to be fine, otherwise ground reflection effects are not well described

% get propagation freq response using ray-tracing (ART)
emission_angle_panam = get_emission_angle(input); % get emission angle from PANAM, to compare with emission angles of the ART
receiver = [ (input{1}.xobs) (input{1}.yobs) (input{1}.zobs) ]; % receiver position
OUT_rayTracing = get_propagation( flight_profile, receiver, nfft, time_PANAM_auralization, emission_angle_panam, show_propagation, tag_auralization) ;

% apply propagation

% consider ground reflection
if isfield( input_file, 'consider_ground_reflection' )
    considerGroundReflection = str2double ( input_file.consider_ground_reflection );   % boolean : 0= only direct path; 1 = direct path + 1st order reflection
else
    considerGroundReflection = 1; % default value
end

% tag_source = 'engineSignal';
% engineSignal = apply_propagation(auralizedEngineSignal, OUT_rayTracing, binaural_signal, show, tag_auralization, tag_source, considerGroundReflection ); % <- binaural truncated to zero
% 
% tag_source = 'airframeSignal';
% airframeSignal = apply_propagation(auralizedAirframeSignal, OUT_rayTracing, binaural_signal, show, tag_auralization, tag_source, considerGroundReflection ); % <- binaural truncated to zero

tag_source = 'overallSignal';
% overallSignal = apply_propagation_FIR1( auralizedOverallSignal, OUT_rayTracing.TF, show, tag_auralization, tag_source, considerGroundReflection );
overallSignal = apply_propagation(auralizedOverallSignal, OUT_rayTracing, binaural_signal, show, tag_auralization, tag_source, considerGroundReflection );

% save final auralized signal (mono signals)
% OutputAuralization.engineSignal = engineSignal.outputSignal;
% OutputAuralization.airframeSignal = airframeSignal.outputSignal;
OutputAuralization.overallSignal = overallSignal.outputSignal;

% attenuation factor (changes dBFS of the written .wav file)
if isfield( input_file, 'AttenuationdB' )
    AttenuationdB = str2double ( input_file.attenuation_db );  % get <AttenuationdB> from <input_file>
else
    AttenuationdB = 0; % default value
end

% save .wav
% fileTag = '_engineSignal';
% save_wav( engineSignal.outputSignal, fs, AttenuationdB, fileTag, tag_auralization );
% 
% fileTag = '_airframeSignal';
% save_wav( airframeSignal.outputSignal, fs, AttenuationdB, fileTag, tag_auralization );

fileTag = '_overallSignal';
save_wav( overallSignal.outputSignal, fs, AttenuationdB, fileTag, tag_auralization );

% binaural signal

if binaural_signal == 1

    % save final auralized signal (stereo signals)
    % OutputAuralization.engineSignal_binaural = engineSignal.outputSignal_binaural;
    % OutputAuralization.airframeSignal_binaural = airframeSignal.outputSignal_binaural;
    OutputAuralization.overallSignal_binaural = overallSignal.outputSignal_binaural;

    % save .wav
    % fileTag = '_engineSignal_binaural';
    % save_wav( engineSignal.outputSignal_binaural, fs, AttenuationdB, fileTag, tag_auralization );
    % 
    % fileTag = '_airframeSignal_binaural';
    % save_wav( airframeSignal.outputSignal_binaural, fs, AttenuationdB, fileTag, tag_auralization );

    fileTag = '_overallSignal_binaural';
    save_wav( overallSignal.outputSignal_binaural, fs, AttenuationdB, fileTag, tag_auralization );

end

%% end

% fprintf('\n- Total auralization (%s-based) time (including processes above):\t%f sec\n', input_type, toc);
fprintf('*--------------------------------------------------------------------------*\n');

end