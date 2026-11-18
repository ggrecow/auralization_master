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
% DOPPLER EFFECT (<doppler_method> in the .ini file, see <get_doppler_method>):
%
%   input (DEFAULT) : the input signal already contains the Doppler effect (sound emission predictions).
%       The FIR filters (both paths, with their relative delay) are applied block by
%       block. The time steps of the predictions are placed at the receiver at constant intervals (dt),
%       i.e. the timing at the receiver follows the timing of the source.
%
%   delay_line : the input signal is synthesized in emission time WITHOUT Doppler
%       effect. For each sound path (direct, reflected):
%       1) a variable delay line reads the signal at the emission time t_e that arrives
%          at each reception time t_r = t_e + tau(t_e), with the propagation time tau of
%          the eigenray (<get_doppler_mapping>, <apply_doppler_delay_line>). This gives
%          the Doppler shift and the compression/stretching of the time axis, and the
%          convective amplification if <convective_amplification> = 1 (DEFAULT 0)
%       2) the FIR filters of the path (without propagation delay, see <get_FIR>) are
%          applied block by block in reception time; each block uses the filter of the
%          time step whose sound is heard at the start of the block
%       3) both paths are summed. Their relative delay comes only from the delay lines.
%       Each time step is thus placed at the time it reaches the receiver, i.e. at uneven
%       intervals. The output starts at the first reception time common to both paths
%       (shifted to zero), and its length differs from the input signal by the change of
%       the propagation time over the flyover.
%       With coherence loss, the decorrelation filters belong to the time steps (epochs
%       of 0.1 s in emission time), so in reception time they change every 0.1*D s.
%
% Author: Gil Felix Greco, Braunschweig  11.03.2024
% 28.03.2025, reformulated to process direct and reflected raypaths
% individually generating a minimum phase FIR (direct path) and a linear phase FIR (reflected path).
% 01.04.2025, included binaural rendering capabilities. The HRTFs are applied individually to the 
% direct and reflected sound paths.
% 16.11.2026, Doppler effect by variable delay lines (<doppler_method> = delay_line)
%
%   AI disclosure (16.11.2026 changes): Generative AI assistance (Anthropic Claude) was used
%   to develop this code under the author's direction. All codes were verified by the author.
% 18.11.2026, optional turbulence-induced amplitude fluctuations applied to the receiver
% signal (mono and binaural), see <apply_amplitude_fluctuations>
%
%   AI disclosure (18.11.2026 changes): Generative AI assistance (Anthropic Claude) was used
%   to develop this code under the author's direction. All codes were verified by the author.
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

%% Doppler effect by variable delay lines (optional, see header)
% the mapping between emission and reception time of each path is computed from the
% propagation times of the time steps of the input data, before they are split for the coherence loss

doppler = get_doppler_method();
useDelayLine = strcmp( doppler.method, 'delay_line' );

if useDelayLine
    nSteps = size( inputRayTracing.propagation_time, 1 );
    t_e = (0:nSteps-1).' * dt_panam; % emission times of the time steps (see <getAuralizationTime>)

    map_direct = get_doppler_mapping( t_e, inputRayTracing.propagation_time(:,1) );
    t_start = map_direct.t_r(1);
    t_end   = map_direct.t_r(end);

    if considerGroundReflection == 1
        map_reflected = get_doppler_mapping( t_e, inputRayTracing.propagation_time(:,2) );
        t_start = max( t_start, map_reflected.t_r(1) );  % first reception time common to both paths
        t_end   = min( t_end,   map_reflected.t_r(end) );
    end

    nOut  = floor( (t_end - t_start)*fs + 1e-6 ) + 1;
    t_out = t_start + (0:nOut-1).' / fs;                 % reception times of the output samples
end

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

if useDelayLine

    % direct path: delay line, then FIR filters in reception time
    nBlocks = floor( nOut/BlockLen );
    y_direct = apply_doppler_delay_line( inputSignal, fs, map_direct, t_out, doppler.convective_amplification );
    idx_direct = il_block_index( map_direct, t_out(1), BlockLen, fs, dt_block, size(FIR.impulseResponse_direct,2), nBlocks );

    outputSignal.outputSignal = overlapp_add_convolution( y_direct, BlockLen, FIR.impulseResponse_direct, 'mono', idx_direct );
    if binaural_signal == 1
        outputSignal.outputSignal_binaural = overlapp_add_convolution( y_direct, BlockLen, FIR.impulseResponse_binaural_direct, 'stereo', idx_direct );
    end

    % reflected path
    if considerGroundReflection == 1
        y_reflected = apply_doppler_delay_line( inputSignal, fs, map_reflected, t_out, doppler.convective_amplification );
        idx_reflected = il_block_index( map_reflected, t_out(1), BlockLen, fs, dt_block, size(FIR.impulseResponse_reflected,2), nBlocks );

        outputSignal.outputSignal = outputSignal.outputSignal + ...
            overlapp_add_convolution( y_reflected, BlockLen, FIR.impulseResponse_reflected, 'mono', idx_reflected );
        if binaural_signal == 1
            outputSignal.outputSignal_binaural = outputSignal.outputSignal_binaural + ...
                overlapp_add_convolution( y_reflected, BlockLen, FIR.impulseResponse_binaural_reflected, 'stereo', idx_reflected );
        end
    end

    % information (printed once, for the overall signal)
    if strcmp( tag_source, 'overallSignal' )
        il_print_doppler_info( map_direct, t_out, numel(inputSignal), fs, doppler, FIR.residualDelay_direct );
        if considerGroundReflection == 1
            fprintf('  reflected path: Doppler factor D = %.3f ... %.3f, propagation-time smoothing rms = %.2g us\n', ...
                min(map_reflected.D), max(map_reflected.D), 1e6*map_reflected.tau_smoothing_rms );
        end
    end

else

    outputSignal.outputSignal = overlapp_add_convolution( inputSignal, BlockLen, FIR.impulseResponse, 'mono' );

    if binaural_signal == 1
        outputSignal.outputSignal_binaural = overlapp_add_convolution( inputSignal, BlockLen, FIR.impulseResponse_binaural, 'stereo' );
    end

end

%% turbulence-induced amplitude fluctuations (optional), see <get_propagation>
% Applied to the receiver signal, i.e. to the sum of the direct and (decorrelated)
% ground-reflected sound, with one random process for all frequencies and both ears
% (<apply_amplitude_fluctuations>). The same seed is used for the engine, airframe and
% overall signals, so that their fluctuations are identical. The DF/FF equalized binaural
% signals are derived from the binaural signal afterwards and contain the same fluctuations.
% Time of each row of sigma at the receiver (middle of each time step, plus the latency of
% the centred FIR filters):
%   input      : constant intervals dt
%   delay_line : arrival time of the direct path, t_e + tau, relative to the output start

if isfield( inputRayTracing, 'amplitudeFluctuations' ) && inputRayTracing.amplitudeFluctuations.consider == 1

    AF      = inputRayTracing.amplitudeFluctuations;
    latency = floor( size( FIR.impulseResponse, 1 )/2 ) / fs;    % IR peak at sample N/2+1
    t_mid   = ( (0:size(AF.sigma,1)-1).' + 0.5 ) * dt_panam;     % middle of each time step (emission)
    if useDelayLine
        t_ctrl = interp1( map_direct.t_e, map_direct.t_r, t_mid, 'linear', 'extrap' ) - t_out(1) + latency;
    else
        t_ctrl = t_mid + latency;
    end
    optsAF = struct( 'seed', AF.seed );

    outputSignal.outputSignal = apply_amplitude_fluctuations( outputSignal.outputSignal, fs, ...
        t_ctrl, AF.f, AF.sigma, optsAF );

    if binaural_signal == 1
        outputSignal.outputSignal_binaural = apply_amplitude_fluctuations( outputSignal.outputSignal_binaural, fs, ...
            t_ctrl, AF.f, AF.sigma, optsAF );
    end

    if strcmp( tag_source, 'overallSignal' )
        fprintf('\n- Amplitude fluctuations applied to the receiver signal(s) (amplitude_fluctuation_seed = %d)\n', AF.seed);
        if show == true
            PLOT_amplitude_fluctuations( AF, t_ctrl, tag_auralization );
        end
    end
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

    % transfer functions with vs. without turbulence-induced coherence loss.
    % The decorrelation filters are identical for the engine, airframe and overall
    % signals (same seed), so they are plotted only once (overall signal)
    if isfield( FIR, 'TF_coherenceLoss' ) && strcmp( tag_source, 'overallSignal' )
        [~, idx_overhead] = min( inputRayTracing.propagation_time(:,1) ); % shortest direct path
        PLOT_coherence_loss_TF( FIR.TF_coherenceLoss, inputRayTracing.coherenceLoss.w, fs, dt_block, ...
            idx_overhead, tag_auralization );
    end
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

%% local function : il_block_index
function idx = il_block_index(map, t0, BlockLen, fs, dt_block, nFilters, nBlocks)
% Filter (column of the FIR matrix) used for each convolution block in reception time:
% the time step whose sound is heard at the start of the block (emission time
% t_e(t_r) at the block start), consistent with <doppler_method> = input, where block k
% of the emission-time signal uses the filter of time step k.
t_block = t0 + (0:nBlocks-1).' * BlockLen/fs;     % reception time at the block starts
te_block = ppval( map.pp_te, t_block );            % emission time heard there
idx = floor( te_block/dt_block + 1e-6 ) + 1;
idx = min( max( idx, 1 ), nFilters );
end

%% local function : il_print_doppler_info
function il_print_doppler_info(map, t_out, nIn, fs, doppler, residualDelay)
fprintf('\n- Doppler effect by variable delay lines (doppler_method = delay_line):\n');
fprintf('  direct path: Doppler factor D = %.3f ... %.3f (f_received = f_emitted/D), propagation-time smoothing rms = %.2g us\n', ...
    min(map.D), max(map.D), 1e6*map.tau_smoothing_rms );
if doppler.convective_amplification == 1
    fprintf('  convective amplification ON: -40*log10(D) = %+.1f ... %+.1f dB\n', ...
        -40*log10(max(map.D)), -40*log10(min(map.D)) );
else
    fprintf('  convective amplification OFF (convective_amplification = 0)\n');
end
fprintf('  duration: %.2f s emitted by the source -> %.2f s at the receiver (first arrival at %.3f s, set to zero)\n', ...
    (nIn-1)/fs, (numel(t_out)-1)/fs, t_out(1) );
maxResidual = max( abs( residualDelay ) );
fprintf('  FIR filters: max. offset of the direct-path peak from the centre = %d sample(s)\n', maxResidual );
if maxResidual > 2
    warning( 'auralization:dopplerResidualDelay', ...
        ['After removing the propagation delay, the direct-path IR peaks up to %d samples away ' ...
         'from the centre. The propagation times and the phase of the transfer functions may not ' ...
         'match; please check the results.'], maxResidual );
end
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