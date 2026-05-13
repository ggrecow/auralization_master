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
%  If binaural files are to be generated, then measured HRTFs
%  from FABIAN database (DOI: 10.14279/depositonce-5718.5)
%  are used (zero-degree HATO only). The aircraft trajectory is already indirectly accounted for
%  by the incidence angle directions, which were previously computed in the <get_propagation>
%  function, and are given as inputs here. Please refer to <AKhrirInterpolation>
%  for the angle convention used by FABIAN (head-centered spherical
%  coordinates).
%
% The procedure to introduce the HRTFs is conducted in frequency domain.
% For that, impulse responses of direct and reflected rays are first alligned
% to the same time delay (centered), and then the relative
% time-of-arrival differences between direct and reflected eigenrays are applied
% to the reflected IRs. This avoids problems related with the periodicity of the
% fourier transform, which leads to a circshift effect for IRs with long time delays,
% and with the overlap and add convolution later on. After that, the impulse
% responses are transformed into freq domain, zero-padded padded as required to perform circular
% convolution, and multiplied to the HRTFs. The resulting
% frequency responses are then transformed into time domain using
% IFFT, resulting in binaural impulse responses for the direct and reflected sound paths individually.
%
% Head orientation is defined by changing the azimuth angles, by the
% <head_orientation> variable read in the input file.  By default, the
% front of the head always points towards the source,
% i.e. positive x-axis -> head_orientation = 0 (degrees), and then angles increase counterclockwise
%  . For ex.,
%   head_orientation = 90 -> source from right ear -> left ear
%   head_orientation = 180 -> source from back  -> to front
%   head_orientation = 270 -> source from left ear -> right ear
%
% Overall, the procedure here is the following, being Hatm the transfer functions containing the atmospheric effects:
%   1) input Hatm is single-sided complex-value spectrum --> get double-sided complex-valued spectrum
%   2) transform double-sided spectrum signal to time domain
%   3) shift impulse response to make it causal (in case of zero-phase FIRs)
%   4) truncate FIR to nTaps. These are FIR filters containing effects
%   described by Hatm
%   5) if both rays are to be considered, then FIR = impulseResponse_direct
%   + impulseResponse_reflected, otherwise FIR = impulseResponse_direct
%
%   If binaural
%   1.b) azimuth and elevation angles are used to get the HRTFs
%   2.b) HRTFs are provided as impulse responses (HRIRs). These are
%   zero-padded to have the same size of the double-sided complex-valued spectrum
%   3.b) zero padded HRIRs are transformed into frequency domain and
%   multiplied with the Hatm in frequency domain
%   4.b) results from step 3.b are transformed into time domain to get
%   impulseResponse_binaural containing both Hatm and HRTFs effects
%   5.b) if both rays are to be considered, then FIR_binaural = impulseResponse_binaural_direct
%   + impulseResponse_binaural_reflected, otherwise FIR_binaural =
%   impulseResponse_direct_binaural
%
%   INPUTS
%
%   inputRayTracing : struct
%   contains the OUT data from ray tracing simulations.
%   type <help get_propagation> for more details
%
%   fs : scalar
%   sampling frquency, in hertz
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
%   caseTag : string
%           flag to sinalize if we are dealing with 'direct' or 'reflected', so that we can save
%           plots with correct name
%
%   output
%   OUT : struct
%       contain following fields
%
%      OUT.impulseResponse : FIR filter truncated to nTaps
%       if  considerGroundReflection=0, OUT.impulseResponse = impulseResponse_direct
%       if  considerGroundReflection=1, OUT.impulseResponse = impulseResponse_direct + impulseResponse_reflected
%
%       OUT.impulseResponse_binaural :  FIR filter truncated to nTaps, for
%       the left and right ears. These FIR contains both the efffects of Hatm and HRTFs together
%       if  considerGroundReflection=0, OUT.impulseResponse_binaural = impulseResponse_direct_binaural
%       if  considerGroundReflection=1, OUT.impulseResponse_binaural = impulseResponse_direct_binaural + impulseResponse_reflected_binaural
%
%       The used HRIRs are provided so that we can plot (or analyse them) them if we need
%       OUT.HRIR_direct = HRIR_direct;
%       OUT.HRIR_reflected= HRIR_reflected;
%
% Gil Felix Greco, Braunschweig 09.04.2025
% Gil Felix Greco, Braunschweig 13.05.2025 - updated IR generation procedure to achieve CORRECT binaural auralizations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global input_file

if binaural_signal == 1
    if isfield( input_file, 'head_orientation' )
        head_orientation = str2double ( input_file.head_orientation );  % get <head_orientation> from <input_file>
    else
        head_orientation = 0; % default value
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

    relativeDelay = idx_r - idx_d;

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

    % direct path
    % incidence angles (spherical coordinates) of the incoming sound wave
    azimuth_direct = input.spherical_angles_HRTF.direct_path(:,1);
    elevation_direct = input.spherical_angles_HRTF.direct_path(:,2);

    % get HRIRs
    [HRIR_direct.leftEar, HRIR_direct.rightEar] = AKhrirInterpolation(azimuth_direct+head_orientation, elevation_direct, HATO, 'measured_sh');
    nSamples_HRIR = size( HRIR_direct.leftEar, 1 );

    % we need to create a zero-padded atmospheric impulse response with size
    % L+M-1 to perform circular convolution with head-related impulse
    % resposes
    atm_IR_direct_padded = [hTOA_direct; zeros(nSamples_HRIR-1, numTimeSteps)];

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

        % we need to create a zero-padded atmospheric impulse response with size
        % L+M-1 to perform circular convolution with head-related impulse
        % resposes
        atm_IR_reflected_padded = [hTOA_reflected; zeros(nSamples_HRIR-1, numTimeSteps)];

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

end