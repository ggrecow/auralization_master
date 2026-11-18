function OUT = get_FIR( input, fs, considerGroundReflection, binaural_signal, tag_auralization )
% function OUT = get_FIR( input, fs, considerGroundReflection, binaural_signal, tag_auralization )
%
% transform input transfer functions into FIR filters. The transfer functions (Hatm) come
% from the ray tracing simulations, and incorporate the following propagation effects
%
%   direct_path : geometrical spreading, atmospheric absorption, phase shift due to propagation time
%
%   reflected_path : geometrical spreading, atmospheric absorption, phase
%   shift due to propagation time, ground reflection
%
% Overall, the procedure is the following (for each time block):
%   1) the single-sided complex-valued Hatm of the direct (and reflected) path is
%      converted to a double-sided spectrum and transformed to time domain (IFFT)
%   2) the direct IR is centered (peak moved to the middle of the IR). This avoids
%      problems related to the periodicity of the Fourier transform, which leads to a
%      circular-shift (wrap-around) effect for IRs with long time delays, and with the
%      overlap-add convolution later on
%   3) if the ground reflection is considered, the reflected IR is centered as well, and
%      the relative time-of-arrival difference between direct and reflected eigenrays
%      (computed on a wrapped circular domain) is applied to the reflected IR only, as a
%      non-circular delay
%   4) OUT.impulseResponse = centered direct IR (+ delayed reflected IR)
%
% With <doppler_method> = delay_line in the .ini file (see <get_doppler_method>), the
% propagation time is applied to each sound path by a variable delay line before the
% FIR filters (Doppler effect, see <apply_propagation>). The FIR filters must then NOT
% contain the propagation delay, and steps 2) to 4) are replaced by:
%   2') the propagation delay of each eigenray is removed from its transfer function,
%       Hatm*exp(-j*2*pi*f*tau) (tau from <get_propagation>; ART defines the delay as
%       exp(+j*2*pi*f*tau), see <AtmosphericPropagation>.PropagationDelayFilter), and the
%       result is conjugated, which converts it to the MATLAB fft sign convention, in
%       which a delay is exp(-j*2*pi*f*tau). The magnitude is unchanged; the conversion
%       keeps the phase of the reflection factor consistent with the real (positive)
%       delay of the reflected path applied by the delay line, so that the sum of both
%       paths has the magnitude of the combined transfer function from ART
%   3') the IR of each path is moved to the centre of the IR (fixed shift of N/2
%       samples, the same for both paths), so both paths have the same latency
%   4') OUT.impulseResponse_direct and OUT.impulseResponse_reflected, applied
%       separately to the delayed signal of each path in <apply_propagation>
%
% If binaural signals are to be generated, measured HRTFs from the FABIAN database
% (DOI: 10.14279/depositonce-5718.5) are used (zero-degree HATO only). The aircraft
% trajectory is indirectly accounted for by the incidence angles of the direct and
% reflected eigenrays, previously computed in <get_propagation> and given as inputs here.
% Please refer to <AKhrirInterpolation> for the angle convention used by FABIAN
% (head-centered spherical coordinates). For each path:
%   1.b) HRIRs are interpolated for the incidence angles (spherical-harmonics
%        representation of the measured HRIRs, 44.1 kHz)
%   2.b) HRIRs are resampled from 44.1 kHz to <fs> (frequency response preserved)
%   3.b) optional high-frequency regularization of the HRIRs (internal flag
%        <hrtf_regularization>, disabled by default; not combinable with
%        <generate_equalized_binaural> = 1 in the .ini file)
%   4.b) atmospheric IRs (from step 2/3) and HRIRs are zero-padded to L+M-1 samples, so
%        that their product in frequency domain corresponds to a linear convolution,
%        and transformed back to time domain
%   5.b) OUT.impulseResponse_binaural = direct binaural IR (+ reflected binaural IR),
%        for the left and right ears
%
% The binaural rendering always uses the measured (blocked ear-canal) HRTFs. The
% diffuse-field (DF) and free-field (FF) equalized versions are obtained afterwards
% from the rendered signal (<generate_equalized_binaural> in the .ini file, see
% <MASTER_auralization_EngineAirframe> and <convert_binaural_equalization>).
%
% Head orientation (<head_orientation> in the .ini file, in degrees, DEFAULT 0, read by
% <get_head_orientation>) is a fixed rotation (yaw) of the head in the horizontal (x-y)
% plane of the PANAM coordinate system. It is the same for all receivers and does not
% follow the aircraft. The HRIRs are interpolated at the azimuth
%   (azimuth of the sound arriving at the receiver) + head_orientation,
% where the arrival azimuth from <get_propagation> is 0 deg = +x and 90 deg = +y (FABIAN
% convention for a head facing +x: 0 = front, 90 = left ear). Increasing values therefore
% turn the head CLOCKWISE seen from above (head facing direction = -head_orientation):
%   head_orientation = 0   -> head faces +x (DEFAULT)
%   head_orientation = 90  -> head faces -y
%   head_orientation = 180 -> head faces -x
%   head_orientation = 270 -> head faces +y
% For a flight towards +x (PANAM approach and departure) and a receiver below the flight
% track (y = 0), the listener hears the aircraft
%   head_orientation = 0   -> from back -> overhead -> front
%   head_orientation = 90  -> from right ear -> left ear
%   head_orientation = 180 -> from front -> overhead -> back (faces the approaching aircraft)
%   head_orientation = 270 -> from left ear -> right ear
% For a flight towards -x, front/back and left/right are swapped. For a lateral receiver,
% the aircraft passes on one side of the head; e.g. head_orientation = 0 and receiver at
% y > 0: back -> right ear -> front. To face the flight track, use 90 (receiver at y > 0)
% or 270 (receiver at y < 0). The arrival and HRTF azimuths of the direct path at the first
% and last time step are printed for each receiver (<print_head_orientation>).
%
%   INPUTS
%
%   input : struct
%   contains the OUT data from ray tracing simulations.
%   type <help get_propagation> for more details
%
%   fs : scalar
%   sampling frequency, in hertz
%
%   considerGroundReflection : boolean
%   contain info about whether to consider ground reflection or not
%   0 = only direct path; 1 = direct path + 1st order reflection
%
%   binaural_signal : logical (boolean)
%   parameter indicating whether binaural signals must be rendered
%
%   tag_auralization : string
%           name to save the plots
%
%   OUTPUTS
%
%   OUT : struct
%       contain following fields
%
%      OUT.impulseResponse : time-varying FIR filters [nSamples x nBlocks]
%       if  considerGroundReflection=0, OUT.impulseResponse = centered direct IR
%       if  considerGroundReflection=1, OUT.impulseResponse = centered direct IR + delayed reflected IR
%
%      OUT.impulseResponse_binaural : time-varying FIR filters for the left and right
%       ears (.leftEar, .rightEar), containing both the effects of Hatm and HRTFs
%       if  considerGroundReflection=0, direct binaural IR
%       if  considerGroundReflection=1, direct binaural IR + reflected binaural IR
%
%      If binaural, the used HRIRs (after the processing steps 1.b-3.b) are also
%      provided, so that they can be plotted or analysed:
%       OUT.HRIR_direct, OUT.HRIR_reflected
%
%      With <doppler_method> = delay_line, the filters of each path are given separately
%      (without propagation delay, centred at sample N/2+1):
%       OUT.impulseResponse_direct, OUT.impulseResponse_reflected [nSamples x nBlocks]
%       OUT.impulseResponse_binaural_direct, OUT.impulseResponse_binaural_reflected
%           (.leftEar, .rightEar), if binaural
%       OUT.residualDelay_direct : [1 x nBlocks] position of the peak of the direct IR
%           relative to the centre, in samples (check: close to 0 if the propagation
%           delay was removed correctly)
%      OUT.impulseResponse and OUT.impulseResponse_binaural (direct + reflected,
%      without relative delay) are then only used for the plots.
%
%      If the turbulence-induced coherence loss is applied, the transfer functions
%      (direct + reflected path, single-sided) without and with the decorrelation
%      filters are also provided, for visualization (see <PLOT_coherence_loss_TF>):
%       OUT.TF_coherenceLoss.coherent, OUT.TF_coherenceLoss.decorrelated
%
% Gil Felix Greco, Braunschweig 09.04.2025
% Gil Felix Greco, Braunschweig 13.05.2026 - updated IR generation procedure to achieve CORRECT binaural auralizations
% Gil Felix Greco, Braunschweig 25.09.2026 - HRTF high-frequency regularization made optional
%   (internal flag <hrtf_regularization>, disabled by default)
% Gil Felix Greco, Braunschweig 25.09.2026 - HRIRs resampled from 44.1 kHz (FABIAN) to <fs>
% Gil Felix Greco, Braunschweig 28.09.2026 - diffuse-field equalization using the FABIAN CTF
%   (<hrtf_df_equalization> in the .ini file, disabled by default)
% Gil Felix Greco, Braunschweig 28.10.2026 - transfer functions with/without coherence loss
%   provided as output (OUT.TF_coherenceLoss) for visualization
% Gil Felix Greco, Braunschweig 11.11.2026 - HRTF equalization removed from the rendering
%   (<hrtf_df_equalization> no longer used): DF and FF equalized binaural signals are now
%   derived from the rendered signal (<generate_equalized_binaural> in the .ini file).
%   Resampled HRIRs scaled by fs_hrir/fs (frequency response preserved)
% Gil Felix Greco, Braunschweig 14.11.2026 - description of <head_orientation> corrected
%   (behaviour unchanged); <head_orientation> read and validated by <get_head_orientation>
% Gil Felix Greco, Braunschweig 16.11.2026 - FIR filters per sound path without
%   propagation delay for <doppler_method> = delay_line
%
%   AI disclosure (11.11.2026, 14.11.2026 and 16.11.2026 changes): Generative AI assistance
%   (Anthropic Claude) was used to develop this code under the author's direction. All
%   codes were verified by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global input_file

if binaural_signal == 1
    head_orientation = get_head_orientation(); % <head_orientation> from <input_file>, DEFAULT 0

    % sampling frequency of the FABIAN HRTF database
    fs_hrir_db = 44100;

    % ---------------------------------------------------------------------
    % OPTIONAL HRTF high-frequency magnitude regularization (see <regularize_hrir>)
    %
    % hrtf_regularization = 0 -> measured HRTFs are used unmodified (DEFAULT)
    % hrtf_regularization = 1 -> HRTF magnitude is flattened (to unity) between
    %                            f1_reg and f2_reg [Hz], phase (ITD) is kept
    %
    % NOTE: this is a perceptual option, not part of the physical rendering.
    % It removes high-frequency spectral cues (pinna notches, high-frequency
    % ILD) relevant for the perception of elevation and front/back direction.
    % ---------------------------------------------------------------------
    hrtf_regularization = 0;
    f1_reg = 4000; % Hz, start of regularization band
    f2_reg = 8000; % Hz, end of regularization band

    if hrtf_regularization == 1
        fprintf('\nHRTF regularization enabled: magnitude flattened between %g Hz and %g Hz\n', f1_reg, f2_reg);
    end

    % The equalized binaural signals (DF, FF; <generate_equalized_binaural> in the .ini)
    % are derived from the rendered signal with filters designed from the UNREGULARIZED
    % HRTFs, so they cannot be combined with <hrtf_regularization>.
    if hrtf_regularization == 1 && isfield( input_file, 'generate_equalized_binaural' ) ...
            && str2double( input_file.generate_equalized_binaural ) == 1
        error( 'auralization:hrtfOptionsConflict', ...
            ['<generate_equalized_binaural> = 1 and <hrtf_regularization> cannot be used together ' ...
            '(the equalization filters were derived from the unregularized HRTFs). Set ' ...
            '<hrtf_regularization> to 0 in get_FIR.m or <generate_equalized_binaural> to 0 in the .ini file.'] );
    end
end

%% Turbulence-induced coherence loss (optional), see <get_propagation>
% The reflected path is filtered with random-phase decorrelation filters
% (Pieren & Lincke 2022). Not applied to flush-mounted receivers.

applyCoherenceLoss = considerGroundReflection == 1 ...
    && isfield( input, 'coherenceLoss' ) && input.coherenceLoss.consider == 1 ...
    && ~( isfield( input, 'flushMountedReceiver' ) && input.flushMountedReceiver );

%% Doppler method (see <get_doppler_method>): with delay lines, the FIR filters of each
% path are designed without the propagation delay

doppler = get_doppler_method();
useDelayLine = strcmp( doppler.method, 'delay_line' );

%% Get atmospheric transfer function and their impulse responses
% direct path

transferFunction =  input.TF;

% get vector dimensions
numTimeSteps = size(transferFunction, 1);
numFreqBins_single_sided = size(transferFunction{1,1}.freq,1);
numFreqBins_double_sided = 2 * (numFreqBins_single_sided - 1);

% Direct path

% declare variable
TF_direct = zeros( numFreqBins_single_sided, numTimeSteps );

for i = 1:numTimeSteps
    TF_direct(:,i) = transferFunction{i,1}.freq(:,1); % [nBins x nTimes]
end

% create double-sided spectra and take its IFFT

% declare variable
TF_direct_double_sided = zeros( numFreqBins_double_sided, numTimeSteps );

for i = 1:numTimeSteps % loop over time steps
    TF_direct_double_sided(:,i) =  [TF_direct(:,i); conj(flipud(TF_direct(2:end-1,i)))];
end

impulseResponse_direct  = ifft(TF_direct_double_sided, 'symmetric');

%% Get atmospheric transfer function and their impulse responses
% Reflected path

if considerGroundReflection == 1

    % declare variables
    TF_reflected = zeros( numFreqBins_single_sided, numTimeSteps );
    TF_combined = zeros( numFreqBins_single_sided, numTimeSteps );

    for i = 1:numTimeSteps
        TF_reflected(:,i) = transferFunction{i,1}.freq(:,2); % [nBins x nTimes]
        TF_combined(:,i) = transferFunction{i,1}.freq(:,3); % only for comparison
    end

    % create double-sided spectra and take its IFFT

    % declare variables
    TF_reflected_double_sided = zeros( numFreqBins_double_sided, numTimeSteps );
    TF_combined_double_sided = zeros( numFreqBins_double_sided, numTimeSteps );

    for i = 1:numTimeSteps % loop over time steps
        TF_reflected_double_sided(:,i) =  [TF_reflected(:,i); conj(flipud(TF_reflected(2:end-1,i)))];
        TF_combined_double_sided(:,i) =  [TF_combined(:,i); conj(flipud(TF_combined(2:end-1,i)))];
    end

    impulseResponse_reflected  = ifft(TF_reflected_double_sided, 'symmetric');
    impulseResponse_combined  = ifft(TF_combined_double_sided, 'symmetric');

    % decorrelated reflected IRs (turbulence-induced coherence loss)
    if applyCoherenceLoss

        % one random phase per filter epoch (0.1 s, see <apply_propagation>):
        % all time steps of the same epoch use the same substream, hence the
        % same random phase (beta still follows the w of each time step).
        % Seeded here, so engine, airframe and overall signals use identical filters.
        if isfield( input.coherenceLoss, 'epoch' )
            epoch = input.coherenceLoss.epoch;
        else
            epoch = (0:numTimeSteps-1).';
        end
        rngStream = RandStream( 'mrg32k3a', 'Seed', input.coherenceLoss.seed );

        impulseResponse_reflected_decorrelated = zeros( size(impulseResponse_reflected) );
        TF_reflected_decorrelated = zeros( size(TF_reflected) ); % kept for <PLOT_coherence_loss_TF>
        for i = 1:numTimeSteps
            rngStream.Substream = epoch(i) + 1;
            B = get_decorrelation_filter( input.coherenceLoss.w(i), numFreqBins_double_sided, fs, rngStream );
            TF_decorrelated = TF_reflected(:,i) .* B;
            TF_reflected_decorrelated(:,i) = TF_decorrelated;
            impulseResponse_reflected_decorrelated(:,i) = ...
                ifft( [TF_decorrelated; conj(flipud(TF_decorrelated(2:end-1)))], 'symmetric' );
        end

    end

end

if useDelayLine

    %% FIR filters per path WITHOUT propagation delay (Doppler effect by delay lines)
    % the delay of each eigenray is removed from its transfer function (ART convention:
    % delay = exp(+j*2*pi*f*tau)) and the result is conjugated (MATLAB convention: delay =
    % exp(-j*2*pi*f*tau)). The zero-delay IRs are moved to the centre (same shift for
    % both paths, i.e. same latency). See the header.

    freq = linspace( 0, fs/2, numFreqBins_single_sided ).';  % as in <get_propagation>
    N = numFreqBins_double_sided;
    centerShift = floor(N/2);                                % peak at sample N/2+1

    tau_direct = input.propagation_time(:,1).';             % [1 x nTimeSteps]
    hTOA_direct = circshift( il_ifft_double_sided( ...
        conj( TF_direct .* exp( -1j*2*pi*freq*tau_direct ) ) ), centerShift );

    % check: the direct IR must peak at the centre
    [~, idxPeak] = max( abs( hTOA_direct ), [], 1 );
    OUT.residualDelay_direct = idxPeak - (centerShift + 1);

    if considerGroundReflection == 1
        tau_reflected = input.propagation_time(:,2).';
        if applyCoherenceLoss
            TF_reflected_used = TF_reflected_decorrelated;
        else
            TF_reflected_used = TF_reflected;
        end
        hTOA_reflected = circshift( il_ifft_double_sided( ...
            conj( TF_reflected_used .* exp( -1j*2*pi*freq*tau_reflected ) ) ), centerShift );
    end

else

%% allign impulse responses of direct path and
% keep relative time-delays of reflected paths from direct paths

% center IRs
[hTOA_direct, ~, idx_d] = il_center_time_varying_ir(impulseResponse_direct);

if considerGroundReflection == 1
    [impulseResponse_reflected_centered, shifts_r, idx_r] = il_center_time_varying_ir(impulseResponse_reflected);

    % with coherence loss: the peak position (and thus the relative delay) is
    % taken from the CLEAN reflected IR; the decorrelated IR gets the same
    % circular shift. The peak of a decorrelated IR is not reliable.
    if applyCoherenceLoss
        for k = 1:numTimeSteps
            impulseResponse_reflected_centered(:,k) = ...
                circshift( impulseResponse_reflected_decorrelated(:,k), shifts_r(k) );
        end
    end

    % compute delay on a wrapped circular domain
    N = size(impulseResponse_direct,1);
    relativeDelay = round( mod(idx_r - idx_d + floor(N/2), N) - floor(N/2) );
    % relativeDelay = idx_r - idx_d;

    % apply relative delays to the IRs of reflected paths
    for k = 1:size(impulseResponse_reflected_centered,2)

        hTOA_reflected(:,k) = delay_ir( ...
            impulseResponse_reflected_centered(:,k), ...
            relativeDelay(k));

    end

else
end

end % if useDelayLine

%% get HRTFs if binaural and apply them to atmospheric transfer functions (in freq domain)
% to get impulseResponses_binaural containing the propagation effect + HRTFs

if binaural_signal == 1

    % head-above-torso orientation, degrees
    % (only one angle is of interest, but more are available in the FABIAN database)
    HATO = 0;

    fs_hrir = fs_hrir_db;

    % direct path
    % incidence angles (spherical coordinates) of the incoming sound wave, i.e. arrival
    % direction at the receiver (0 = +x, 90 = +y); + head_orientation gives the head-related
    % azimuth (head turned clockwise by head_orientation, see header)
    azimuth_direct = input.spherical_angles_HRTF.direct_path(:,1);
    elevation_direct = input.spherical_angles_HRTF.direct_path(:,2);

    % get HRIRs
    [HRIR_direct.leftEar, HRIR_direct.rightEar] = AKhrirInterpolation(azimuth_direct+head_orientation, elevation_direct, HATO, 'measured_sh');

    % resample HRIRs from the FABIAN sampling frequency (44.1 kHz) to <fs>
    HRIR_direct.leftEar  = il_resample_hrir( HRIR_direct.leftEar,  fs, fs_hrir );
    HRIR_direct.rightEar = il_resample_hrir( HRIR_direct.rightEar, fs, fs_hrir );
    nSamples_HRIR = size( HRIR_direct.leftEar, 1 );

    % we need to create a zero-padded atmospheric impulse response with size
    % L+M-1 to perform circular convolution with head-related impulse
    % resposes
    atm_IR_direct_padded = [hTOA_direct; zeros(nSamples_HRIR-1, numTimeSteps)];

    % optional high-frequency regularization of the HRIRs
    if hrtf_regularization == 1
        HRIR_direct.leftEar  = regularize_hrir( HRIR_direct.leftEar, fs, f1_reg, f2_reg );
        HRIR_direct.rightEar = regularize_hrir( HRIR_direct.rightEar, fs, f1_reg, f2_reg );
    end

    % zero-pad HRIRs so they have L+M-1 samples
    HRIR_direct_leftEar_padded = [HRIR_direct.leftEar; zeros(numFreqBins_double_sided - 1, numTimeSteps)];
    HRIR_direct_rightEar_padded = [HRIR_direct.rightEar; zeros(numFreqBins_double_sided - 1, numTimeSteps)];

    % get binaural impulse responses (i.e. atmospheric effects + HRTFs)
    impulseResponse_direct_binaural_leftEar = ifft( fft(atm_IR_direct_padded) .* fft( HRIR_direct_leftEar_padded ), 'symmetric' );
    impulseResponse_direct_binaural_rightEar = ifft( fft(atm_IR_direct_padded) .* fft( HRIR_direct_rightEar_padded ), 'symmetric' );

    % reflected path
    if considerGroundReflection == 1

        % incidence angle (spherical coordinates) of the incoming sound wave
        azimuth_reflected = input.spherical_angles_HRTF.reflected_path(:,1);
        elevation_reflected = input.spherical_angles_HRTF.reflected_path(:,2);

        % get HRIRs for reflected path
        [HRIR_reflected.leftEar, HRIR_reflected.rightEar] = AKhrirInterpolation(azimuth_reflected+head_orientation, elevation_reflected, HATO, 'measured_sh');

        % resample HRIRs from the FABIAN sampling frequency (44.1 kHz) to <fs>
        HRIR_reflected.leftEar  = il_resample_hrir( HRIR_reflected.leftEar,  fs, fs_hrir );
        HRIR_reflected.rightEar = il_resample_hrir( HRIR_reflected.rightEar, fs, fs_hrir );

        % we need to create a zero-padded atmospheric impulse response with size
        % L+M-1 to perform circular convolution with head-related impulse
        % resposes
        atm_IR_reflected_padded = [hTOA_reflected; zeros(nSamples_HRIR-1, numTimeSteps)];

        % optional high-frequency regularization of the HRIRs
        if hrtf_regularization == 1
            HRIR_reflected.leftEar  = regularize_hrir( HRIR_reflected.leftEar, fs, f1_reg, f2_reg );
            HRIR_reflected.rightEar = regularize_hrir( HRIR_reflected.rightEar, fs, f1_reg, f2_reg );
        end

        % zero-pad HRIRs so they have L+M-1 samples
        HRIR_reflected_leftEar_padded = [HRIR_reflected.leftEar; zeros(numFreqBins_double_sided - 1, numTimeSteps)];
        HRIR_reflected_rightEar_padded = [HRIR_reflected.rightEar; zeros(numFreqBins_double_sided - 1, numTimeSteps)];

        % get binaural impulse responses (i.e. atmospheric effects + HRTFs)
        impulseResponse_reflected_binaural_leftEar = ifft( fft(atm_IR_reflected_padded) .* fft( HRIR_reflected_leftEar_padded ), 'symmetric' );
        impulseResponse_reflected_binaural_rightEar = ifft( fft(atm_IR_reflected_padded) .* fft( HRIR_reflected_rightEar_padded), 'symmetric' );

    end
end

%% get output variables

if useDelayLine
    % filters of each path, applied separately after the delay lines (<apply_propagation>)
    OUT.impulseResponse_direct = hTOA_direct;
    if binaural_signal == 1
        OUT.impulseResponse_binaural_direct.leftEar  = impulseResponse_direct_binaural_leftEar;
        OUT.impulseResponse_binaural_direct.rightEar = impulseResponse_direct_binaural_rightEar;
    end
    if considerGroundReflection == 1
        OUT.impulseResponse_reflected = hTOA_reflected;
        if binaural_signal == 1
            OUT.impulseResponse_binaural_reflected.leftEar  = impulseResponse_reflected_binaural_leftEar;
            OUT.impulseResponse_binaural_reflected.rightEar = impulseResponse_reflected_binaural_rightEar;
        end
    end
end

if considerGroundReflection == 1

    % OUT.impulseResponse = impulseResponse_combined ;
    OUT.impulseResponse = ( hTOA_direct + hTOA_reflected );

    if binaural_signal == 1

        % variables used for visualization only (see <PLOT_FIR> and <PLOT_FIR_spectrogram> below)
        OUT.impulseResponse_binaural.leftEar =  ( impulseResponse_direct_binaural_leftEar + impulseResponse_reflected_binaural_leftEar );
        OUT.impulseResponse_binaural.rightEar =  ( impulseResponse_direct_binaural_rightEar + impulseResponse_reflected_binaural_rightEar );

        % HRIRs from FABIAN - used for visualizytion only
        OUT.HRIR_direct = HRIR_direct;
        OUT.HRIR_reflected= HRIR_reflected;
    end

    % transfer functions without and with turbulence-induced coherence loss
    % (single-sided, [nBins x nTimeSteps]) - used for visualization only (see <PLOT_coherence_loss_TF>).
    % Both are built as direct + reflected, so they differ ONLY by the decorrelation
    % filters actually applied in the auralization (same seed and 0.1 s epochs)
    if applyCoherenceLoss
        OUT.TF_coherenceLoss.coherent     = TF_direct + TF_reflected;
        OUT.TF_coherenceLoss.decorrelated = TF_direct + TF_reflected_decorrelated;
    end

else % only direct path required

    OUT.impulseResponse = hTOA_direct;

    if binaural_signal == 1

        OUT.impulseResponse_binaural.leftEar = ( impulseResponse_direct_binaural_leftEar );
        OUT.impulseResponse_binaural.rightEar = ( impulseResponse_direct_binaural_rightEar );

        % HRIRs from FABIAN - used for visualizytion only
        OUT.HRIR_direct = HRIR_direct;
    end

end


%% plot freq response of FIR filter

tBlock = round(numTimeSteps/2); % a time bin to plot the freq response of the FIR filter

if useDelayLine

    % filters of each path (without propagation delay): magnitude compared with the
    % transfer function of the same path
    PLOT_FIR( OUT.impulseResponse_direct, TF_direct, tBlock, fs, [tag_auralization '_direct'] );
    PLOT_FIR_spectrogram( OUT.impulseResponse_direct, TF_direct, fs, [tag_auralization '_direct'] );

    if considerGroundReflection == 1
        PLOT_FIR( OUT.impulseResponse_reflected, TF_reflected, tBlock, fs, [tag_auralization '_reflected'] );
        PLOT_FIR_spectrogram( OUT.impulseResponse_reflected, TF_reflected, fs, [tag_auralization '_reflected'] );
    end

    if binaural_signal == 1
        PLOT_FIR( OUT.impulseResponse_binaural_direct.leftEar, TF_direct, tBlock, fs, [tag_auralization '_direct_binaural_left'] );
        PLOT_FIR_spectrogram( OUT.impulseResponse_binaural_direct.leftEar, TF_direct, fs, [tag_auralization '_direct_binaural_left'] );
    end

elseif considerGroundReflection == 1

    % plot FIR for a fixed (single) source/receiver combination
    PLOT_FIR( OUT.impulseResponse, TF_combined, 52, fs, [tag_auralization '_combined'] );

    % plot FIR spectrogram
    PLOT_FIR_spectrogram(OUT.impulseResponse, TF_combined, fs, [tag_auralization '_combined'] );

    if binaural_signal == 1

        % plot FIR for a fixed (single) source/receiver combination
        PLOT_FIR( OUT.impulseResponse_binaural.leftEar, TF_combined, tBlock, fs, [tag_auralization '_combined_binaural_left'] );

        % plot FIR spectrogram
        PLOT_FIR_spectrogram(OUT.impulseResponse_binaural.leftEar, TF_combined, fs, [tag_auralization '_combined_binaural_left'] );

    end

else

    % plot FIR for a fixed (single) source/receiver combination
    PLOT_FIR( OUT.impulseResponse, TF_direct, tBlock, fs, [tag_auralization '_direct'] );

    % plot FIR spectrogram
    PLOT_FIR_spectrogram(OUT.impulseResponse, TF_direct, fs, [tag_auralization '_direct'] );

    if binaural_signal == 1

        % plot FIR for a fixed (single) source/receiver combination
        PLOT_FIR( OUT.impulseResponse_binaural.leftEar, TF_direct, tBlock, fs, [tag_auralization '_direct_binaural_left'] );

        % plot FIR spectrogram
        PLOT_FIR_spectrogram(OUT.impulseResponse_binaural.leftEar, TF_direct, fs, [tag_auralization '_direct_binaural_left'] );

    end

end

%% inline function : il_ifft_double_sided(TF)

    function IR_out = il_ifft_double_sided(TF_in)
        % single-sided spectra [nBins x nTimes] -> real IRs [2*(nBins-1) x nTimes]
        IR_out = ifft( [TF_in; conj(flipud(TF_in(2:end-1,:)))], 'symmetric' );
    end % end of function <il_ifft_double_sided>

%% inline function : il_resample_hrir(HRIR, fs, fs_hrir)

    function HRIR_out = il_resample_hrir(HRIR_in, fs, fs_hrir)
        % resample HRIRs [Nsamples x Nblocks] from <fs_hrir> to <fs>.
        % The FABIAN database is sampled at 44.1 kHz; without resampling, the
        % HRIRs would be played back at the framework sampling frequency,
        % shifting all spectral features by fs/fs_hrir and scaling all ITDs
        % by fs_hrir/fs.
        %
        % <resample> preserves the AMPLITUDE of a signal, i.e. the resampled IR has
        % fs/fs_hrir times more samples of (about) the same values, so its frequency
        % response would be scaled by fs/fs_hrir (+0.74 dB at 48 kHz). For an impulse
        % response, the result is therefore multiplied by fs_hrir/fs, which keeps the
        % frequency response (e.g. 0 dB at low frequencies) unchanged.
        if fs == fs_hrir
            HRIR_out = HRIR_in;
        else
            HRIR_out = resample( HRIR_in, fs, fs_hrir ) * ( fs_hrir/fs ); % operates along columns
        end
    end % end of function <il_resample_hrir>

%% inline function : il_center_time_varying_ir(IR)

    function [IR_centered, shifts, peakIndices] = il_center_time_varying_ir(IR)

        % CENTER_TIME_VARYING_IR
        %
        % Centers each time-varying impulse response independently
        % using circular shift based on maximum energy.
        %
        % INPUT
        %   IR : [Nsamples x Nblocks]
        %
        % OUTPUT
        %   IR_centered : centered IR matrix
        %   shifts      : applied shifts per block
        %   peakIndices : peak indices per block

        % dimensions
        [N, Nblocks] = size(IR);

        % allocate outputs
        IR_centered = zeros(size(IR));
        shifts      = zeros(1, Nblocks);
        peakIndices = zeros(1, Nblocks);

        % center sample index
        centerIndex = floor(N/2) + 1;

        % loop over time-varying IRs
        for k = 1:Nblocks

            % current IR
            h = IR(:,k);

            % energy
            energy = abs(h).^2;

            % maximum-energy index
            [~, peakIdx] = max(energy);

            % required shift
            shift = centerIndex - peakIdx;

            % store
            shifts(k) = shift;
            peakIndices(k) = peakIdx;

            % centered IR
            IR_centered(:,k) = circshift(h, shift);

        end

    end % end of function <il_center_time_varying_ir>

%% inline function : delay_ir(h, delaySamples)

    function h_delayed = delay_ir(h, delaySamples)

        % DELAY_IR
        %
        % Applies a PHYSICAL (non-circular) delay to an impulse response
        % in the time domain.
        %
        % Unlike MATLAB's circshift(), this function DOES NOT wrap samples
        % around the vector boundaries. Samples shifted beyond the vector
        % length are discarded and replaced with zeros.
        %
        % This behavior is physically meaningful for acoustic propagation,
        % since delayed energy should not reappear at the beginning of the IR.
        %
        %
        % -------------------------------------------------------------------------
        % INPUTS
        % -------------------------------------------------------------------------
        %
        % h : [N x 1] vector
        %
        %     Input impulse response in the TIME DOMAIN.
        %
        %     The function assumes:
        %       - h is causal or approximately causal
        %       - rows correspond to time samples
        %
        %
        % delaySamples : scalar integer
        %
        %     Desired delay in samples.
        %
        %     Positive values:
        %         Delay the IR toward later times
        %         (shift RIGHT)
        %
        %     Negative values:
        %         Advance the IR toward earlier times
        %         (shift LEFT)
        %
        %     Examples:
        %
        %         delaySamples = 20
        %
        %             h[n] -> h[n-20]
        %
        %             Adds 20 zeros at the beginning.
        %
        %
        %         delaySamples = -10
        %
        %             h[n] -> h[n+10]
        %
        %             Removes first 10 samples and pads zeros at end.
        %
        %
        % -------------------------------------------------------------------------
        % OUTPUT
        % -------------------------------------------------------------------------
        %
        % h_delayed : [N x 1] vector
        %
        %     Delayed impulse response with SAME LENGTH as input.
        %
        %     The vector length is preserved by:
        %
        %       - truncating samples shifted outside the vector
        %       - zero-padding empty regions
        %
        %
        % -------------------------------------------------------------------------
        % IMPORTANT NOTES
        % -------------------------------------------------------------------------
        %
        % 1) INTEGER-SAMPLE DELAYS ONLY
        %
        %     This function only supports integer delays.
        %
        %     If sub-sample delays are required, then:
        %
        %         - interpolation
        %         - fractional delay filters
        %         - or frequency-domain phase shifts
        %
        %     should be used instead.
        %
        % -------------------------------------------------------------------------
        %
        % 4) SIGN CONVENTION
        %
        %     Positive delay:
        %
        %         shifts energy toward larger sample indices
        %
        %         --> later arrival
        %
        %     Negative delay:
        %
        %         shifts energy toward smaller sample indices
        %
        %         --> earlier arrival
        %
        % -------------------------------------------------------------------------

        h = h(:);

        N = length(h);

        if delaySamples >= 0

            h_delayed = [zeros(delaySamples,1); h];

            h_delayed = h_delayed(1:N);

        else

            delaySamples = abs(delaySamples);

            h_delayed = h(delaySamples+1:end);

            h_delayed(end+1:N) = 0;

        end

    end

%%

    function HRIR_out = regularize_hrir(HRIR_in, fs, f1, f2)
        % REGULARIZE_HRIR
        %
        % Smoothly regularizes the HIGH-FREQUENCY MAGNITUDE RESPONSE of
        % Head-Related Impulse Responses (HRIRs) while preserving:
        %
        %   - interaural time differences (ITD)
        %   - phase structure
        %   - group delay
        %   - low-frequency interaural level differences (ILD)
        %
        % -------------------------------------------------------------------------
        % MOTIVATION
        % -------------------------------------------------------------------------
        %
        % Raw free-field HRTFs often contain very sharp spectral peaks and notches
        % at high frequencies (> ~6 kHz), mainly caused by pinna resonances.
        %
        % For long-distance aircraft noise auralization, these detailed spectral
        % structures may become overly prominent and can:
        %
        %   - amplify weak high-frequency numerical noise
        %   - exaggerate coloration
        %   - produce metallic or artificial timbre
        %   - reduce perceptual realism
        %
        % This function progressively FLATTENS the HRTF MAGNITUDE RESPONSE above
        % a chosen frequency range while preserving the ORIGINAL PHASE RESPONSE.
        %
        % As a consequence:
        %
        %   - binaural timing cues remain intact
        %   - localization is largely preserved
        %   - excessive spectral coloration is reduced
        %
        % -------------------------------------------------------------------------
        % PRINCIPLE
        % -------------------------------------------------------------------------
        %
        % The complex HRTF is decomposed as:
        %
        %       H(f) = |H(f)| exp(j*phi(f))
        %
        % Only the magnitude |H(f)| is modified.
        %
        % The phase phi(f) is preserved unchanged.
        %
        % The modified magnitude is:
        %
        %       |H_reg(f)| =
        %
        %           (1-alpha(f)) |H(f)| + alpha(f)
        %
        % where:
        %
        %       alpha(f) = 0      below f1
        %       alpha(f) = 1      above f2
        %
        % and transitions smoothly between f1 and f2.
        %
        % Therefore:
        %
        %   below f1:
        %       original HRTF preserved
        %
        %   above f2:
        %       magnitude becomes flat (= 1)
        %
        % -------------------------------------------------------------------------
        % INPUTS
        % -------------------------------------------------------------------------
        %
        % HRIR_in : [Nsamples x Nblocks]
        %
        %       Matrix containing HRIRs.
        %
        %       Each column corresponds to one HRIR.
        %
        %       Example:
        %
        %           HRIR_leftEar
        %           HRIR_rightEar
        %
        %
        % fs : scalar
        %
        %       Sampling frequency [Hz]
        %
        %
        % f1 : scalar
        %
        %       Start frequency of regularization [Hz]
        %
        %       Below f1:
        %           HRTF remains unchanged.
        %
        %
        % f2 : scalar
        %
        %       End frequency of regularization [Hz]
        %
        %       Above f2:
        %           magnitude response becomes flat.
        %
        %
        % -------------------------------------------------------------------------
        % OUTPUT
        % -------------------------------------------------------------------------
        %
        % HRIR_out : [Nsamples x Nblocks]
        %
        %       Regularized HRIRs.
        %
        % -------------------------------------------------------------------------

        %% dimensions

        [N, Nblocks] = size(HRIR_in);

        %% allocate output

        HRIR_out = zeros(size(HRIR_in));

        %% frequency vector

        f = (0:N-1)'/N * fs;

        %% build smooth transition function alpha(f)

        alpha = zeros(size(f));

        % fully regularized above f2
        alpha(f > f2) = 1;

        % smooth cosine transition between f1 and f2
        idx_transition = (f >= f1) & (f <= f2);

        alpha(idx_transition) = ...
            0.5 * ( ...
            1 - cos( ...
            pi * (f(idx_transition)-f1)/(f2-f1) ...
            ));

        %% process each HRIR independently

        for k = 1:Nblocks

            %% FFT of HRIR

            H = fft(HRIR_in(:,k));

            %% magnitude and phase

            magnitude = abs(H);

            phase = angle(H);

            %% regularize ONLY magnitude

            magnitude_reg = ...
                (1 - alpha) .* magnitude + alpha;

            %% reconstruct complex HRTF

            H_reg = magnitude_reg .* exp(1j * phase);

            %% enforce Hermitian symmetry
            %
            % This guarantees a real-valued HRIR after IFFT.
            %
            % Only required for numerical robustness.

            if mod(N,2) == 0

                % even length FFT

                H_reg(N/2+2:end) = ...
                    conj(flipud(H_reg(2:N/2)));

            else

                % odd length FFT

                H_reg((N+3)/2:end) = ...
                    conj(flipud(H_reg(2:(N+1)/2)));

            end

            %% back to time domain

            HRIR_out(:,k) = real(ifft(H_reg));

        end

    end

end