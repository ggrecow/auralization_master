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
% Gil Felix Greco, Braunschweig 14.11.2026 - head orientation check printed for each
%   receiver (<print_head_orientation>)
% Gil Felix Greco, Braunschweig 16.11.2026 - Doppler effect by variable delay lines
%   (<doppler_method> = delay_line, see <get_doppler_method>); Doppler factor from the
%   propagation times compared with the PANAM straight-line estimate
%
%   AI disclosure (14.11.2026 and 16.11.2026 changes): Generative AI assistance (Anthropic
%   Claude) was used to develop this code under the author's direction. All codes were
%   verified by the author.
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

% DF and FF equalized versions of the binaural signal (only if binaural_signal = 1).
% Read (and validated) here, before the time-consuming rendering. See the example .ini files.
if isfield( input_file, 'generate_equalized_binaural' )
    generate_equalized_binaural = str2double( input_file.generate_equalized_binaural );
    if ~( generate_equalized_binaural == 0 || generate_equalized_binaural == 1 )
        error( 'auralization:invalidGenerateEqualizedBinaural', ...
            '<generate_equalized_binaural> in the .ini file must be 0 or 1 (got %s).', ...
            input_file.generate_equalized_binaural );
    end
else
    generate_equalized_binaural = 0; % default value
end

% keys of previous (unreleased) versions are no longer used
if isfield( input_file, 'hrtf_df_equalization' )
    fprintf(['\n- NOTE: <hrtf_df_equalization> in the .ini file is no longer used and is ignored. ' ...
             'The binaural signal is always rendered with the measured HRTFs; use ' ...
             '<generate_equalized_binaural> = 1 to also obtain DF and FF equalized versions.\n']);
end

% Doppler method (see <get_doppler_method>), read and validated before the rendering.
% It is a decision of the user and must match the PANAM input data (with or without
% Doppler effect), which cannot be checked here.
doppler = get_doppler_method();
if strcmp( doppler.method, 'delay_line' )
    if ~strcmp( input_type, 'emission' )
        error( 'auralization:dopplerInputType', ...
            '<doppler_method> = delay_line can only be used with emission data (input_type = emission).' );
    end
    if doppler.convective_amplification == 1
        txt_ca = 'ON';
    else
        txt_ca = 'OFF';
    end
    fprintf(['\n- Doppler effect: variable delay lines (doppler_method = delay_line). The PANAM ' ...
             'emission data must NOT contain the Doppler effect. Convective amplification: %s\n'], txt_ca );
else
    fprintf('\n- Doppler effect: included in the PANAM emission data (doppler_method = input)\n');
    if isfield( input_file, 'convective_amplification' )
        fprintf('  NOTE: <convective_amplification> is only used with doppler_method = delay_line and is ignored.\n');
    end
end

% time resolution from PANAM. Assumes it is constant (source noise prediction dt from PANAM is usually .5 seconds and cte) - ACHTUNG: this is the case only for sound source   
global dt_panam
dt_panam =  input{2,1}.source_time - input{1,1}.source_time;  

% check the sampling frequency (before the time-consuming rendering)
il_check_sampling_frequency( fs, dt_panam, binaural_signal );

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
show_propagation = show; % plot propagation results

% FFT length used for freq vector of atmospheric transfer
% function and as nTaps of the FIR filter of atmospheric transfer function
% nfft = round(fs*dt_panam);
% nfft = str2double( input_file.nfft );
nfft = 2^14; % to yield df=3 Hz (if fs=48 kHz). needs to be fine, otherwise ground reflection effects are not well described

% get propagation freq response using ray-tracing (ART)
emission_angle_panam = get_emission_angle(input); % get emission angle from PANAM, to compare with emission angles of the ART
receiver = [ (input{1}.xobs) (input{1}.yobs) (input{1}.zobs) ]; % receiver position

% A flush-mounted (ground-plane) receiver and a head-and-torso simulator are mutually
% exclusive: at zobs = 0 there is no head, and the HRTFs would be meaningless.
% The mono signal is the physically correct output in this case, and it already
% contains the +6 dB from p_tot = p_dir*(1+Q). See issue #4.
if receiver(3) == 0 && binaural_signal == 1
    warning(['Flush-mounted receiver (zobs = 0 m): binaural rendering disabled. ' ...
             'Only the mono (ground-plane) signal is generated. ' ...
             'Use zobs = 1.2-1.5 m for a HATS receiver.']);
    binaural_signal = 0;
end

OUT_rayTracing = get_propagation( flight_profile, receiver, nfft, time_PANAM_auralization, emission_angle_panam, show_propagation, tag_auralization) ;

% Doppler factor from the propagation times vs. PANAM straight-line estimate (check only)
if strcmp( doppler.method, 'delay_line' )
    il_check_doppler_factor( input, OUT_rayTracing, dt_panam, show, tag_auralization );
end

% print the HRTF angles of the direct path at the first and last time step (console and
% log.txt), so the head orientation (<head_orientation> in the .ini file) can be checked
if binaural_signal == 1
    print_head_orientation( OUT_rayTracing.spherical_angles_HRTF, get_head_orientation(), receiver );
end

% apply propagation

% consider ground reflection
if isfield( input_file, 'consider_ground_reflection' )
    considerGroundReflection = str2double ( input_file.consider_ground_reflection );   % boolean : 0= only direct path; 1 = direct path + 1st order reflection
else
    considerGroundReflection = 1; % default value
end

tag_source = 'engineSignal';
engineSignal = apply_propagation(auralizedEngineSignal, OUT_rayTracing, binaural_signal, show, tag_auralization, tag_source, considerGroundReflection ); % <- binaural truncated to zero

tag_source = 'airframeSignal';
airframeSignal = apply_propagation(auralizedAirframeSignal, OUT_rayTracing, binaural_signal, show, tag_auralization, tag_source, considerGroundReflection ); % <- binaural truncated to zero

tag_source = 'overallSignal';
% overallSignal = apply_propagation_FIR1( auralizedOverallSignal, OUT_rayTracing.TF, show, tag_auralization, tag_source, considerGroundReflection );
overallSignal = apply_propagation(auralizedOverallSignal, OUT_rayTracing, binaural_signal, show, tag_auralization, tag_source, considerGroundReflection );

% save final auralized signal (mono signals)
OutputAuralization.engineSignal = engineSignal.outputSignal;
OutputAuralization.airframeSignal = airframeSignal.outputSignal;
OutputAuralization.overallSignal = overallSignal.outputSignal;

% attenuation factor (changes dBFS of the written .wav file)
if isfield( input_file, 'attenuation_db' )
    AttenuationdB = str2double ( input_file.attenuation_db );  % get <attenuation_db> from <input_file>
else
    AttenuationdB = 0; % default value
end

% save .wav
fileTag = '_engineSignal';
save_wav( engineSignal.outputSignal, fs, AttenuationdB, fileTag, tag_auralization );

fileTag = '_airframeSignal';
save_wav( airframeSignal.outputSignal, fs, AttenuationdB, fileTag, tag_auralization );

fileTag = '_overallSignal';
save_wav( overallSignal.outputSignal, fs, AttenuationdB, fileTag, tag_auralization );

% binaural signal

if binaural_signal == 1 % <- binaural truncated to zero for engineSignal and airframeSignal 

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

    % optional: diffuse-field (DF) and free-field (FF) equalized versions of the SAME
    % rendering (see <convert_binaural_equalization>; exact, since both equalizations are
    % fixed linear filters per ear and the block convolution is linear):
    %   *_overallSignal_binaural_DFeq_FS... -> playback over diffuse-field equalized headphones
    %   *_overallSignal_binaural_FFeq_FS... -> computing metrics (LAeq, EPNL, ECMA-418-2 'free-frontal', ...)
    if generate_equalized_binaural == 1

        OutputAuralization.overallSignal_binaural_DFeq = convert_binaural_equalization( ...
            overallSignal.outputSignal_binaural, fs, 'none', 'DF' );
        fileTag = '_overallSignal_binaural_DFeq';
        save_wav( OutputAuralization.overallSignal_binaural_DFeq, fs, AttenuationdB, fileTag, tag_auralization );

        OutputAuralization.overallSignal_binaural_FFeq = convert_binaural_equalization( ...
            overallSignal.outputSignal_binaural, fs, 'none', 'FF' );
        fileTag = '_overallSignal_binaural_FFeq';
        save_wav( OutputAuralization.overallSignal_binaural_FFeq, fs, AttenuationdB, fileTag, tag_auralization );

        fprintf('\n- Diffuse-field (DF) and free-field (FF) equalized binaural signals saved\n');
    end

end

%% end

% fprintf('\n- Total auralization (%s-based) time (including processes above):\t%f sec\n', input_type, toc);
fprintf('*--------------------------------------------------------------------------*\n');

end

%% local function: check sampling frequency
function il_check_doppler_factor( input, OUT_rayTracing, dt_panam, show, tag_auralization )
% Doppler factor D = dt_r/dt_e of the direct path, as used by the delay lines (from the
% propagation times of the eigenrays), compared with the straight-line estimate from the
% PANAM data, D = 1 - M*cos(phix), used for the buzzsaw tones with doppler_method = input
% (see <getTonalInput>). Both should agree closely for a homogeneous atmosphere without
% wind; refraction and wind change D from the propagation times only.

global save_mat_fig

nSteps = size( OUT_rayTracing.propagation_time, 1 );
t_e = (0:nSteps-1).' * dt_panam;
map = get_doppler_mapping( t_e, OUT_rayTracing.propagation_time(:,1) );

hasPanam = all( isfield( input{1}, {'vel', 'sound_speed', 'phix'} ) );
if hasPanam
    D_panam = zeros( nSteps, 1 );
    for i = 1:nSteps
        D_panam(i) = 1 - ( input{i}.vel ./ input{i}.sound_speed ) .* cos( deg2rad( input{i}.phix ) );
    end
    fprintf(['- Doppler factor of the direct path: %.3f ... %.3f (from propagation times), ' ...
             '%.3f ... %.3f (PANAM straight-line estimate); max. difference %.3f\n'], ...
        min(map.D), max(map.D), min(D_panam), max(D_panam), max( abs( map.D - D_panam ) ) );
end

if show == true
    figure('name', 'PROCESSING - Doppler factor (direct path)');
    h = gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

    a = plot( t_e, map.D, 'k' ); hold on;
    leg = { 'from propagation times (delay line)' };
    if hasPanam
        b = plot( t_e, D_panam, 'r--' );
        a = [a b];
        leg{end+1} = 'PANAM: $1 - M\cos\phi_x$';
    end
    grid on;
    xlabel('Emission time, $t_{\mathrm{e}}$ (s)', 'Interpreter', 'Latex');
    ylabel('Doppler factor, $D = \mathrm{d}t_{\mathrm{r}}/\mathrm{d}t_{\mathrm{e}}$ (-)', 'Interpreter', 'Latex');
    legend( a, leg, 'Interpreter', 'Latex', 'Location', 'SE' );
    set(gcf,'color','w');

    if ~isempty( tag_auralization )
        filename = strcat( tag_auralization, '_Doppler_factor' );
        save_pdf = 1; save_png = 0;
        export_figures( filename, save_mat_fig, save_png, save_pdf );
    end
end
end

function il_check_sampling_frequency( fs, dt_panam, binaural_signal )
% Checks <sampling_freq> from the .ini file:
%   - fs >= 44.1 kHz (error otherwise): the HRIRs (FABIAN, 44.1 kHz) and the
%     synthesis would lose high-frequency content, which biases the metrics
%   - fs*dt_panam must be an integer number of samples (error otherwise): the
%     block length is round(fs*dt_panam), so a non-integer value would make the
%     auralized signal drift slowly in time over the flyover
%   - note if fs > 48 kHz: no benefit (binaural content ends at 22.05 kHz and air
%     absorption removes the rest), only longer computation
%   - note if binaural and fs ~= 44.1 kHz: HRIRs are resampled, so the binaural
%     signal contains no energy above 22.05 kHz (the mono signal does)
% Recommended values: 48000 Hz (DEFAULT; native rate of ECMA-418-2 and of most
% audio hardware) or 44100 Hz.

if fs < 44100 || fs ~= round( fs )
    error( 'auralization:invalidSamplingFrequency', ...
        ['<sampling_freq> in the .ini file must be an integer >= 44100 Hz (got %g Hz). ' ...
         'Recommended: 48000 Hz (DEFAULT) or 44100 Hz.'], fs );
end

nBlock = fs*dt_panam;
if abs( nBlock - round( nBlock ) ) > 1e-3   % tolerance for rounding errors in the PANAM time vector
    error( 'auralization:nonIntegerBlockLength', ...
        ['<sampling_freq> = %g Hz and the PANAM time step dt = %g s give a non-integer block ' ...
         'length of %.4f samples, which would make the auralized signal drift in time. ' ...
         'Use a sampling frequency for which fs*dt is an integer (e.g. 48000 or 44100 Hz).'], ...
         fs, dt_panam, nBlock );
end

if fs > 48000
    fprintf( ['\n- NOTE: <sampling_freq> = %g Hz. Rates above 48 kHz bring no benefit (binaural ' ...
              'content ends at 22.05 kHz, air absorption removes the rest) and increase the ' ...
              'computation time.\n'], fs );
end

if binaural_signal == 1 && fs ~= 44100
    fprintf( ['\n- NOTE: <sampling_freq> = %g Hz: the FABIAN HRIRs (44.1 kHz) are resampled. The ' ...
              'binaural signals contain no energy above 22.05 kHz, the mono signal does (up to ' ...
              '%g kHz). Irrelevant for A-weighted levels, EPNL and loudness; consider it when ' ...
              'comparing unweighted levels of mono and binaural signals.\n'], fs, fs/2000 );
end

end
