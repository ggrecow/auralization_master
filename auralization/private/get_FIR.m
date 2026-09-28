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
% If binaural signals are to be generated, measured HRTFs from the FABIAN database
% (DOI: 10.14279/depositonce-5718.5) are used (zero-degree HATO only). The aircraft
% trajectory is indirectly accounted for by the incidence angles of the direct and
% reflected eigenrays, previously computed in <get_propagation> and given as inputs here.
% Please refer to <AKhrirInterpolation> for the angle convention used by FABIAN
% (head-centered spherical coordinates). For each path:
%   1.b) HRIRs are interpolated for the incidence angles (spherical-harmonics
%        representation of the measured HRIRs, 44.1 kHz)
%   2.b) optional diffuse-field equalization (<hrtf_df_equalization> in the .ini file,
%        enabled by default): HRIRs are filtered with the inverted common transfer
%        function (CTF) of the FABIAN database, at 44.1 kHz
%   3.b) HRIRs are resampled from 44.1 kHz to <fs>
%   4.b) optional high-frequency regularization of the HRIRs (internal flag
%        <hrtf_regularization>, disabled by default; not combinable with 2.b)
%   5.b) atmospheric IRs (from step 2/3) and HRIRs are zero-padded to L+M-1 samples, so
%        that their product in frequency domain corresponds to a linear convolution,
%        and transformed back to time domain
%   6.b) OUT.impulseResponse_binaural = direct binaural IR (+ reflected binaural IR),
%        for the left and right ears
%
% Head orientation is defined by changing the azimuth angles, by the
% <head_orientation> variable read in the input file. By default, the front of the
% head points to the positive x-axis (head_orientation = 0 degrees), and angles
% increase counterclockwise. For ex.,
%   head_orientation = 0   -> source from front -> back (DEFAULT)
%   head_orientation = 90  -> source from right ear -> left ear
%   head_orientation = 180 -> source from back -> front
%   head_orientation = 270 -> source from left ear -> right ear
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
%      If binaural, the used HRIRs (after the processing steps 1.b-4.b) are also
%      provided, so that they can be plotted or analysed:
%       OUT.HRIR_direct, OUT.HRIR_reflected
%
% Gil Felix Greco, Braunschweig 09.04.2025
% Gil Felix Greco, Braunschweig 13.05.2026 - updated IR generation procedure to achieve CORRECT binaural auralizations
% Gil Felix Greco, Braunschweig 25.09.2026 - HRTF high-frequency regularization made optional
%   (internal flag <hrtf_regularization>, disabled by default)
% Gil Felix Greco, Braunschweig 25.09.2026 - HRIRs resampled from 44.1 kHz (FABIAN) to <fs>
% Gil Felix Greco, Braunschweig 28.09.2026 - diffuse-field equalization using the FABIAN CTF
%   (<hrtf_df_equalization> in the .ini file, enabled by default)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global input_file

if binaural_signal == 1
    if isfield( input_file, 'head_orientation' )
        head_orientation = str2double ( input_file.head_orientation );  % get <head_orientation> from <input_file>
    else
        head_orientation = 0; % default value
    end

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

    % diffuse-field (DF) equalization of the HRTFs, defined in the .ini file
    % (see <hrtf_df_equalization> in the .ini for a description)
    %
    % hrtf_df_equalization = 1 -> HRIRs are convolved with the inverted, 3rd-octave
    %                             smoothed, minimum-phase common transfer function (CTF)
    %                             of the FABIAN database (DEFAULT). Same processing as
    %                             AKhrirInterpolation(..., 'dir').
    % hrtf_df_equalization = 0 -> measured (blocked ear-canal) HRTFs are used
    %
    % Uses <FABIAN_CTF_measured_inverted_smoothed.mat>, shipped with the framework
    % in third_party/FABIAN_HRTF_DATABASE_v4/ (extracted from the FABIAN database
    % with utilities/convert_FABIAN_CTF_sofa2mat.m).
    % Not to be combined with <hrtf_regularization>.
    if isfield( input_file, 'hrtf_df_equalization' )
        hrtf_df_equalization = str2double( input_file.hrtf_df_equalization );
    else
        hrtf_df_equalization = 0; % default value
    end

    if ~( hrtf_df_equalization == 0 || hrtf_df_equalization == 1 )
        error( 'auralization:invalidDfEqualization', ...
            '<hrtf_df_equalization> in the .ini file must be 0 or 1 (got %s).', ...
            num2str( hrtf_df_equalization ) );
    end

    if hrtf_df_equalization == 1

        if hrtf_regularization == 1
            error( 'auralization:hrtfOptionsConflict', ...
                ['<hrtf_df_equalization> and <hrtf_regularization> cannot be used together ' ...
                '(the CTF was derived from the unregularized HRTFs). Set <hrtf_regularization> ' ...
                'to 0 in get_FIR.m or <hrtf_df_equalization> to 0 in the .ini file.'] );
        end

        ctf_file = 'FABIAN_CTF_measured_inverted_smoothed.mat';
        if ~exist( ctf_file, 'file' )
            error( 'auralization:ctfNotFound', ...
                ['Diffuse-field equalization requires <%s>, which is shipped in ' ...
                'third_party/FABIAN_HRTF_DATABASE_v4/. Make sure this folder is on the MATLAB path.'], ctf_file );
        end

        CTF = load( ctf_file, 'ctf', 'fs' );

        if CTF.fs ~= fs_hrir_db
            error( 'auralization:ctfSamplingFrequency', ...
                'CTF sampling frequency (%g Hz) differs from the HRTF database (%g Hz).', CTF.fs, fs_hrir_db );
        end

        fprintf('\nHRTF diffuse-field equalization enabled (FABIAN inverted CTF)\n');
    else
        fprintf('\nHRTF diffuse-field equalization disabled: measured (blocked ear-canal) HRTFs are used\n');
    end
end

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

end

%% allign impulse responses of direct path and
% keep relative time-delays of reflected paths from direct paths

% center IRs
[hTOA_direct, ~, idx_d] = il_center_time_varying_ir(impulseResponse_direct);

if considerGroundReflection == 1
    [impulseResponse_reflected_centered, ~, idx_r] = il_center_time_varying_ir(impulseResponse_reflected);

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

%% get HRTFs if binaural and apply them to atmospheric transfer functions (in freq domain)
% to get impulseResponses_binaural containing the propagation effect + HRTFs

if binaural_signal == 1

    % head-above-torso orientation, degrees
    % (only one angle is of interest, but more are available in the FABIAN database)
    HATO = 0;

    fs_hrir = fs_hrir_db;

    % direct path
    % incidence angles (spherical coordinates) of the incoming sound wave
    azimuth_direct = input.spherical_angles_HRTF.direct_path(:,1);
    elevation_direct = input.spherical_angles_HRTF.direct_path(:,2);

    % get HRIRs
    [HRIR_direct.leftEar, HRIR_direct.rightEar] = AKhrirInterpolation(azimuth_direct+head_orientation, elevation_direct, HATO, 'measured_sh');

    % optional diffuse-field equalization (at the database sampling frequency,
    % i.e. before resampling; same operation as AKhrirInterpolation(..., 'dir'))
    if hrtf_df_equalization == 1
        HRIR_direct.leftEar  = fftfilt( CTF.ctf, HRIR_direct.leftEar );
        HRIR_direct.rightEar = fftfilt( CTF.ctf, HRIR_direct.rightEar );
    end

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

        % optional diffuse-field equalization (at the database sampling frequency,
        % i.e. before resampling; same operation as AKhrirInterpolation(..., 'dir'))
        if hrtf_df_equalization == 1
            HRIR_reflected.leftEar  = fftfilt( CTF.ctf, HRIR_reflected.leftEar );
            HRIR_reflected.rightEar = fftfilt( CTF.ctf, HRIR_reflected.rightEar );
        end

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

if considerGroundReflection == 1

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

%% inline function : il_resample_hrir(HRIR, fs, fs_hrir)

    function HRIR_out = il_resample_hrir(HRIR_in, fs, fs_hrir)
        % resample HRIRs [Nsamples x Nblocks] from <fs_hrir> to <fs>.
        % The FABIAN database is sampled at 44.1 kHz; without resampling, the
        % HRIRs would be played back at the framework sampling frequency,
        % shifting all spectral features by fs/fs_hrir and scaling all ITDs
        % by fs_hrir/fs.
        if fs == fs_hrir
            HRIR_out = HRIR_in;
        else
            HRIR_out = resample( HRIR_in, fs, fs_hrir ); % operates along columns
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