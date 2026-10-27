function outputSignal = apply_propagation(inputSignal, inputRayTracing, binaural_signal, show, tag_auralization, tag_source, considerGroundReflection )
% function outputSignal = apply_propagation(inputSignal, inputRayTracing, binaural_signal, show, tag_auralization, tag_source, considerGroundReflection )
%
% The auralization of the aircraft noise at a receiver position is conducted here by applying the 
% atmospheric <transferFunction> computed using the fuction <get_propagation> to the <inputSignal> 
% (sinthesized sound source sound file). The transfer functions are
% transformed into FIR filters of a truncated length nTaps and applied to
% the <inputSignal> using uniformelly partitioned circular convolution.
%
% If required, generates stereo (binaural) signals using measured HRTFs 
% from FABIAN database (zero-degree HATO) (DOI:
% 10.14279/depositonce-5718.5) . See <get_FIR.m> function for more details. 
% 
% INPUT:
%       inputSignal : [N x 1] vector
%       Auralized source description, in Pascal 
%
%       inputRayTracing : struct
%       contains the following fields/data from ray tracing simulations.
%       type <help get_propagation> for more details
%     
%       binaural_signal : logical (boolean)
%       parameter indicating whether binaural signals must be rendered 
%
%       show : logical (boolean)
%       optional parameter for figures (results) display
%       'false' (disable, default value) or 'true' (enable)
%
%       tag_source : string
%       provide info whether the signal being considered comes from the airframe, 
%       engine or overall noise emission. Used only on the show/save plots
%
%       considerGroundReflection : boolean
%       contain info about whether to consider ground reflection or not
%       0 = only direct path; 1 = direct path + 1st order reflection 
%
%       OUTPUT : struct
%
%       outputSignal : [Ntime x 1] vector
%       Auralized aircraft noise in the ground, in Pascal (mono signal)
%
%       outputSignal_binaural : [Ntime x 2] vector
%       Auralized aircraft noise in the ground, in Pascal (stereo signal)
%
% Author: Gil Felix Greco, Braunschweig  11.03.2024
% 28.03.2025, reformulated to process direct and reflected raypaths
% individually generating a minimum phase FIR (direct path) and a linear phase FIR (reflected path).
% 01.04.2025, included binaural rendering capabilities. The HRTFs are applied individually to the 
% direct and reflected sound paths.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global fs % sampling frequency used for auralization, defined in the <auralization_master>
global dt_panam % dt from panam, defined in the <auralization_master>

switch tag_source

    case 'engineSignal'
        binaural_signal = 0;

    case 'airframeSignal'
        binaural_signal = 0;
end

if isfield(inputRayTracing, 'flushMountedReceiver') && inputRayTracing.flushMountedReceiver
    binaural_signal = 0;
end

tic;

%% get transform atmospheric transfer function as FIR filters
% and head-related impulse responses (if required)

%% turbulence-induced coherence loss: decorrelation filters redrawn every 0.1 s
% (Pieren & Lincke 2022), independent of dt_panam. Filters can only change at
% block boundaries, so each PANAM step is split onto the common sample grid of
% dt_panam and 0.1 s. Within a PANAM step the deterministic propagation is held
% (as before); only the random decorrelation phase changes.

dt_block = dt_panam;

useCoherenceLoss = considerGroundReflection == 1 ...
    && isfield( inputRayTracing, 'coherenceLoss' ) && inputRayTracing.coherenceLoss.consider == 1 ...
    && ~( isfield( inputRayTracing, 'flushMountedReceiver' ) && inputRayTracing.flushMountedReceiver );

if useCoherenceLoss

    t_update = 0.1;                         % [s] filter update interval
    nSubMax  = 20;                          % cost limit for the sub-block split

    L_block  = round( fs*dt_panam );        % samples per PANAM step
    L_update = round( fs*t_update );        % samples per filter epoch
    L_sub    = gcd( L_block, L_update );    % common grid of both
    nSub     = L_block / L_sub;             % sub-blocks per PANAM step

    if nSub <= nSubMax
        if nSub > 1
            inputRayTracing = il_split_time_steps( inputRayTracing, nSub );
            dt_block = L_sub / fs;
        end
        k = ( 0:numel(inputRayTracing.TF)-1 ).';
        inputRayTracing.coherenceLoss.epoch = floor( k*L_sub / L_update );   % exact 0.1 s epochs
        fprintf('\n- Coherence loss: decorrelation filters redrawn every %.3f s (%d sub-block(s) per PANAM step)\n', ...
            t_update, nSub);
    else
        % no practical common grid: each block uses the epoch in which it starts
        k = ( 0:numel(inputRayTracing.TF)-1 ).';
        inputRayTracing.coherenceLoss.epoch = floor( k*L_block / L_update );
        warning( 'auralization:coherenceUpdate', ...
            ['dt_panam = %.4f s and 0.1 s have no practical common sample grid (fs = %d Hz). ' ...
             'Decorrelation filters are redrawn at the block boundary following each 0.1 s mark.'], ...
            dt_panam, fs );
    end
end

FIR = get_FIR( inputRayTracing, fs, considerGroundReflection, binaural_signal, tag_auralization );

%% apply FIRs to aircraft (emission+Doppler) signal using circular convolution

% define synthesis block length [samples]
BlockLen = round(fs*dt_block);

outputSignal.outputSignal = overlapp_add_convolution( inputSignal, BlockLen, FIR.impulseResponse, 'mono' );

if binaural_signal == 1
    outputSignal.outputSignal_binaural = overlapp_add_convolution( inputSignal, BlockLen, FIR.impulseResponse_binaural, 'stereo' );
end

% %% check plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% plot(outputSignal.outputSignal );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot signals

if show == true
    tag_title =  ['OUTPUT - Spectrogram of auralized emission after propagation - ' tag_source];
    tag_save = ['_' tag_source '_spectrogram_emission_propagated'];
    PLOT_spectrogram(outputSignal.outputSignal, fs, tag_title, tag_auralization, tag_save);
else
end

% sound(outputSignal.outputSignal, 44100, 24 );

if binaural_signal == 1
    if show == true
        % binaural - spectrogram
        tag_title =  ['OUTPUT - Spectrogram of auralized emission after propagation (binaural) - ' tag_source];
        tag_save = ['_' tag_source '_spectrogram_emission_propagated_binaural'];
        PLOT_spectrogram(outputSignal.outputSignal_binaural, fs, tag_title, tag_auralization, tag_save);

        switch tag_source
            case 'overallSignal'

                % binaural - spectrogram of applied HRTFs
                % (reflected HRIRs only exist if consider_ground_reflection = 1)
                if isfield( FIR, 'HRIR_reflected' )
                    HRIR_reflected = FIR.HRIR_reflected;
                else
                    HRIR_reflected = []; % only the direct path is plotted
                end

                tag_title =  ['OUTPUT - Spectrogram of  HRIRs  - ' tag_source];
                tag_save = ['_' tag_source '_spectrogram_HRIR'];
                PLOT_HRIR_spectrogram(FIR.HRIR_direct, HRIR_reflected, dt_block, tag_title, tag_auralization, tag_save);

                %  spectrogram - difference between mono and stereo (binaural) signals
                tag_title =  ['OUTPUT - Spectrogram (binaural - mono) - ' tag_source];
                tag_save = ['_' tag_source '_spectrogram_binaural_mono'];
                PLOT_spectrogram_difference(outputSignal.outputSignal, outputSignal.outputSignal_binaural, fs, tag_title, tag_auralization, tag_save);
        end

    end
end

fprintf('\n- Applying propagation FIR filter/auralization time:\t%f sec\n', toc);

end

%% local function : il_split_time_steps

function in = il_split_time_steps(in, nSub)
% Repeat every PANAM time step nSub times. Deterministic propagation (TFs,
% HRTF angles, propagation times, coherence parameter w) is held within a
% PANAM step; only the random decorrelation phase changes between sub-blocks
% (see <get_FIR>).

in.TF = repelem( in.TF, nSub, 1 );
in.spherical_angles_HRTF.direct_path    = repelem( in.spherical_angles_HRTF.direct_path,    nSub, 1 );
in.spherical_angles_HRTF.reflected_path = repelem( in.spherical_angles_HRTF.reflected_path, nSub, 1 );
in.propagation_time = repelem( in.propagation_time, nSub, 1 );
in.coherenceLoss.w  = repelem( in.coherenceLoss.w,  nSub, 1 );

end
