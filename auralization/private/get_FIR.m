function OUT = get_FIR( input, fs, considerGroundReflection, binaural_signal, tag_auralization )
% function OUT = get_FIR( input, fs, considerGroundReflection, binaural_signal, tag_auralization )
%
% transform input transfer functions into FIR filters. The transfer functions (Hatm) come 
% from ray tracing simulations, and incorporate the following propagation effects 
%
%   direct_path : geometrical spreading, atmospheric absorption, phase shift due to propagation time
%
%   reflected_path : geometrical spreading, atmospheric absorption, phase
%   shift due to propagation time, ground reflection
%
%  If binaural files are to be generated, then measured HRTFs 
%  from FABIAN database (zero-degree HATO) (DOI: 10.14279/depositonce-5718.5) 
%  are used. The aircraft trajectory is already indirectly accounted for  
%  by the incidence angle directions, which were previously computed in the <get_propagation> 
%  function, and are given as inputs here. Please refer to <AKhrirInterpolation>
%  for the angle convention used by FABIAN (head-centered spherical
%  coordinates). The procedure to introduce the HRTFs is conducted in frequency domain.
% The HRTFs are obtained for the direct and reflected rays individually, and then multiplied to 
% respective atmospheric transfer functions in frequency domain. The
% resulting frequency responses are then transformed into time domain using
% IFFT. 
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
% The FIR filters (considering or not HRTFs) contain both direct and
% reflected rays. This is done in order to avoid loosing the phase (or
% relative delays) between the impulse responses of the direct and
% reflected rays. Alligning the relative time-shifts of FIR filters of
% direct and reflected rays has been proven (to me) to be a hassle because
% the phase of the reflected ray is also influenced by freq-dependent
% effects caused by ground reflections, and not only by propagation time
% delay. IF those phase changes are not accounted (somehow) when generating
% the FIR filters, then overlap and add convolution does not incorporate
% all propagation effeects correctly. 
%
% Overall, the procedure here is the following, being Hatm the transfer functions containing the atmospheric effects: 
%   1) input Hatm is single-sided complex-value spectrum --> get double-sided complex-valued spectrum
%   2) transform double-sided spectrum signal to time domain
%   3) shift impulse response to make it causal
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
% OBS: it is chosen to define (truncate) the FIR nTaps so that a lower frequency of 10
% Hz is achieved. This is because low frequencies are important for
% aircraft noise.
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

transferFunction =  input.TF;

% get vector dimensions
numTimeSteps = size(transferFunction, 1);
numFreqBins_single_sided = size(transferFunction{1,1}.freq,1);
numFreqBins_double_sided = 2 * (numFreqBins_single_sided - 1);  

% Direct path
TF_direct = zeros( numFreqBins_single_sided, numTimeSteps ); % declare variable

for i = 1:numTimeSteps
    TF_direct(:,i) = transferFunction{i,1}.freq(:,1); % [nBins x nTimes]
end

% create double-sided spectra and take its IFFT
% TF_direct_double_sided = AKsingle2bothSidedSpectrum(TF_direct);

TF_direct_double_sided = zeros( numFreqBins_double_sided, numTimeSteps ); % declare variable

for i = 1:numTimeSteps % loop over time steps
    TF_direct_double_sided(:,i) =  [TF_direct(:,i); conj(flipud(TF_direct(2:end-1,i)))];
end

impulseResponse_direct  = ifft(TF_direct_double_sided, 'symmetric');

% guarantee causal impulseResponses
impulseResponse_direct = il_makeCausal (impulseResponse_direct );

% Reflected path
if considerGroundReflection == 1

    TF_reflected =zeros( numFreqBins_single_sided, numTimeSteps ); % declare variable
    TF_combined =zeros( numFreqBins_single_sided, numTimeSteps );

    for i = 1:numTimeSteps
        TF_reflected(:,i) = transferFunction{i,1}.freq(:,2); % [nBins x nTimes]
        TF_combined(:,i) = transferFunction{i,1}.freq(:,3); % only for comparison
    end

    % create double-sided spectra and take its IFFT
    % TF_reflected_double_sided = AKsingle2bothSidedSpectrum(TF_reflected);

    TF_reflected_double_sided = zeros( numFreqBins_double_sided, numTimeSteps ); % declare variable
    TF_combined_double_sided = zeros( numFreqBins_double_sided, numTimeSteps ); % declare variable
    for i = 1:numTimeSteps % loop over time steps
        TF_reflected_double_sided(:,i) =  [TF_reflected(:,i); conj(flipud(TF_reflected(2:end-1,i)))];
        TF_combined_double_sided(:,i) =  [TF_combined(:,i); conj(flipud(TF_combined(2:end-1,i)))];
    end

    impulseResponse_reflected  = ifft(TF_reflected_double_sided, 'symmetric');
    impulseResponse_combined  = ifft(TF_combined_double_sided, 'symmetric');

    % guarantee causal impulseResponses
    impulseResponse_reflected = il_makeCausal (impulseResponse_reflected );
    impulseResponse_combined= il_makeCausal (impulseResponse_combined );

end

%% get HRTFs if binaural and apply them to atmospheric transfer functions (in freq domain)
% to get impulseResponses_binaural containing the propagation effect + HRTFs  

if binaural_signal == 1

    % direct path
    % incidence angle (spherical coordinates) of the incoming sound wave (direct path)
    azimuth_direct = input.spherical_angles_HRTF.direct_path(:,1);
    elevation_direct = input.spherical_angles_HRTF.direct_path(:,2);
    HATO = 0; % head-above-torso orientation, degrees

    % get HRIRs - direct path
    [HRIR_direct.l_a, HRIR_direct.r_a] = AKhrirInterpolation(azimuth_direct+head_orientation, elevation_direct, HATO, 'measured_sh');
    nSamples_HRIR = size( HRIR_direct.l_a, 1);

    % zero-pad HRIRs so they have the same size as the atmospheric transfer functions
    HRIR_direct_leftEar_padded = [HRIR_direct.l_a; zeros(numFreqBins_double_sided - nSamples_HRIR, numTimeSteps)];
    HRIR_direct_rightEar_padded = [HRIR_direct.r_a; zeros(numFreqBins_double_sided - nSamples_HRIR, numTimeSteps)];

    % get TF_binaural (i.e. atmospheric effects + HRTFs)
    impulseResponse_direct_binaural_leftEar = ifft( TF_direct_double_sided .* fft( HRIR_direct_leftEar_padded ), 'symmetric' );
    impulseResponse_direct_binaural_rightEar = ifft( TF_direct_double_sided .* fft( HRIR_direct_rightEar_padded), 'symmetric' );

    % guarantee causal impulseResponses
    impulseResponse_direct_binaural_leftEar = il_makeCausal ( impulseResponse_direct_binaural_leftEar );
    impulseResponse_direct_binaural_rightEar = il_makeCausal ( impulseResponse_direct_binaural_rightEar );

    % reflected path
    if considerGroundReflection == 1

        % incidence angle (spherical coordinates) of the incoming sound wave (direct path)
        azimuth_reflected = input.spherical_angles_HRTF.reflected_path(:,1);
        elevation_reflected = input.spherical_angles_HRTF.reflected_path(:,2);

        % get HRIRs - reflected path
        [HRIR_reflected.l_a, HRIR_reflected.r_a] = AKhrirInterpolation(azimuth_reflected+head_orientation, elevation_reflected, HATO, 'measured_sh');

        % zero-pad HRIRs so they have the same size as the atmospheric transfer functions
        HRIR_reflected_leftEar_padded = [HRIR_reflected.l_a; zeros(numFreqBins_double_sided - nSamples_HRIR, numTimeSteps)];
        HRIR_reflected_rightEar_padded = [HRIR_reflected.r_a; zeros(numFreqBins_double_sided - nSamples_HRIR, numTimeSteps)];

        % get TF_binaural (i.e. atmospheric effects + HRTFs)
        impulseResponse_reflected_binaural_leftEar = ifft( TF_reflected_double_sided .* fft( HRIR_reflected_leftEar_padded ), 'symmetric' );
        impulseResponse_reflected_binaural_rightEar = ifft( TF_reflected_double_sided .* fft( HRIR_reflected_rightEar_padded), 'symmetric' );

        % guarantee causal impulseResponses
        impulseResponse_reflected_binaural_leftEar = il_makeCausal (impulseResponse_reflected_binaural_leftEar );
        impulseResponse_reflected_binaural_rightEar = il_makeCausal (impulseResponse_reflected_binaural_rightEar );

    end
end

%% truncate FIR to nTaps 

% define nTaps
min_freq = 10; % in (Hz)
nTaps = 2^14 ;

if considerGroundReflection == 1

    OUT.impulseResponse = il_truncate_FIR ( impulseResponse_combined, nTaps);
    % OUT.impulseResponse = il_truncate_FIR ( impulseResponse_direct + impulseResponse_reflected, nTaps);

    if binaural_signal == 1
        OUT.impulseResponse_binaural.left_ear = il_truncate_FIR ( impulseResponse_direct_binaural_leftEar + impulseResponse_reflected_binaural_leftEar, nTaps);
        OUT.impulseResponse_binaural.right_ear = il_truncate_FIR ( impulseResponse_direct_binaural_rightEar + impulseResponse_reflected_binaural_rightEar, nTaps);
        OUT.HRIR_direct = HRIR_direct;
        OUT.HRIR_reflected= HRIR_reflected;
    end

else

    OUT.impulseResponse = il_truncate_FIR ( impulseResponse_direct, nTaps);

    if binaural_signal == 1
        OUT.impulseResponse_binaural.left_ear = il_truncate_FIR ( impulseResponse_direct_binaural_leftEar, nTaps);
        OUT.impulseResponse_binaural.right_ear = il_truncate_FIR ( impulseResponse_direct_binaural_rightEar, nTaps);
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
        PLOT_FIR( OUT.impulseResponse_binaural.left_ear, TF_combined, tBlock, fs, [tag_auralization '_combined_binaural_left'] );

        % plot FIR spectrogram
        PLOT_FIR_spectrogram(OUT.impulseResponse_binaural.left_ear, TF_combined, fs, [tag_auralization '_combined_binaural_left'] );

    end

else

    % plot FIR for a fixed (single) source/receiver combination
    PLOT_FIR( OUT.impulseResponse, TF_direct, tBlock, fs, [tag_auralization '_direct'] );

    % plot FIR spectrogram
    PLOT_FIR_spectrogram(OUT.impulseResponse, TF_direct, fs, [tag_auralization '_direct'] );

    if binaural_signal == 1

        % plot FIR for a fixed (single) source/receiver combination
        PLOT_FIR( OUT.impulseResponse_binaural.left_ear, TF_direct, tBlock, fs, [tag_auralization '_direct_binaural_left'] );

        % plot FIR spectrogram
        PLOT_FIR_spectrogram(OUT.impulseResponse_binaural.left_ear, TF_direct, fs, [tag_auralization '_direct_binaural_left'] );

    end

end

%% function: make causal impulse response

    function output = il_makeCausal (input )
    % function output = il_makeCausal (input )
    % make impulse reponse causal by shifting M/2 samples (if necessary)
    %
    % More info: https://www.youtube.com/watch?v=R6sKlHHTULM&list=LL&index=12&t=766s
    %
    % input - impulse response [nSamples,:]
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % shift impulse response to make it simmetrical
    Nsamples = size(input,1); % number of samples of the impulse response

     % sometimes it is necessary to mirror the IFFT(FRF), sometimes not. This loop
     % tries to deal with that in a very coarse way. This problem/solution to be further
     % investigated/improved
     for i =1:size(input,2)

         if input(end,i) > 1e-8 % if this arbitrary threshold is achieved in the last bin, then it means we need to mirror the IR
             output(:,i) = [input(Nsamples/2+1:end,i); input(1:Nsamples/2,i)  ];
         else
             output(:,i) = input(:,i);
         end

     end

    end % end of function <il_makeCausal>

%% truncate FIR to nTaps

    function output = il_truncate_FIR ( input, nTaps)
    % function output = il_truncate_FIR ( input, nTaps)
    %
    % truncate impulse response to nTaps
    %
    % input - impulse response [nSamples,:]
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    numCols = size(input,2);

    % zero pad IR to avoid problems when windowing FIR to a desired number of taps
    input_padded = [ zeros(nTaps, numCols); input; zeros(nTaps, numCols)];

    for i = 1:numCols

        [~, idxMax(i)] = max ( input_padded(:,i) );
        output(:,i) = input_padded( idxMax(i) - (nTaps/2) : idxMax(i) + (nTaps/2) - 1, i); % ntaps/2

    end

    end % end of function <il_truncate_FIR>

end