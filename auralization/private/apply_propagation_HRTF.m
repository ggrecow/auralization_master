function outputSignal = apply_propagation_HRTF(inputSignal, inputRayTracing, show, tag_auralization, tag_source, considerGroundReflection )
% function outputSignal = apply_propagation_HRTF(inputSignal, transferFunction, show, tag_auralization, tag_source, considerGroundReflection)
%
% The auralization of the aircraft noise at a receiver position is conducted here by applying the 
% atmospheric <transferFunction> computed using the fuction <get_propagation> to the <inputSignal> 
% (sinthesized sound source sound file). The transfer functions are
% transformed into FIR filter of a truncated length nTaps and applied to
% the <inputSignal> using uniformelly partitioned circular convolution.
%
% INPUT:
%       inputSignal : [N x 1] vector
%       Auralized source description, in Pascal 
%
%       inputRayTracing : struct
%       contains the following fields/data from ray tracing simulations.
%       type <help get_propagation> for more details
%     
%       nTaps : scalar
%       nfft or number of freq bins considered in the FIR filter. Dictates freq discretization, and 
%       thus also the min freq that can be considered by the FIR filters     
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
% OUTPUT:
%       outputSignal : [N x 1] vector
%       Auralized aircraft noise in the ground, in Pascal
%
% Author: Gil Felix Greco, Braunschweig  11.03.2024
% 28.03.2025, reformulated to process direct and reflected raypaths
% individually. Include capability of generating binaural signals.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global fs % sampling frequency used for auralization, defined in the <auralization_master>
global dt_panam % dt from panam, defined in the <auralization_master>

tic;
%% Inputs

transferFunction = inputRayTracing.TF;
propagation_time = inputRayTracing.propagation_time;
spherical_angles_HRTF = inputRayTracing.spherical_angles_HRTF;

ray_delay = propagation_time(:,2) - propagation_time(:,1); % delay between reflected and direct ray paths, in seconds
ray_delay_samples = round ( ( propagation_time(:,2) .* fs) - ( propagation_time(:,1) .* fs) ); % delay between reflected and direct ray paths, in seconds

%% Get atmospheric transfer function 

 % declare variables to pre allocate memory
nTF = size(transferFunction,1);
TF_direct = zeros( size(transferFunction,1), size(transferFunction{1,1}.freq,1) );
TF_reflection = zeros( size(transferFunction,1), size(transferFunction{1,1}.freq,1) );
% TF_combined = zeros( size(transferFunction,1), size(transferFunction{1,1}.freq,1) );

for i = 1:nTF
    TF_direct(i,:) = transferFunction{i,1}.freq(:,1); % size [nTimes x nBins]
    TF_reflection(i,:) = transferFunction{i,1}.freq(:,2);
    % TF_combined(i,:) = transferFunction{i,1}.freq(:,3);
end

% freq vector from FR - single sided spectrum - nBins = (round(fs*dt_panam)+1)/2 (see <get_propagation> function)
% freq = transferFunction{1,1}.freqVector(:,1);

% number of blocks (i.e. time instances)
% nTime = size(transferFunction,1);  

clear nTF

%% transform atmospheric transfer function into FIR filter

% direct sound path - lets work only with nTimes in the columns
transferFunction_singleSideSpectrum_direct = transpose( TF_direct ); % single sided spectrum - nBins = (round(fs*dt_panam)+1)/2 (see <get_propagation> function)

caseTag = 'direct';
FIR_direct = get_FIR_from_atm_TF( transferFunction_singleSideSpectrum_direct, ray_delay_samples, fs, tag_auralization, caseTag );

if considerGroundReflection == 1

   % reflected sound path - lets work only with nTimes in the columns
    transferFunction_singleSideSpectrum_reflected = transpose( TF_reflection ); % single sided spectrum - nBins = (round(fs*dt_panam)+1)/2 (see <get_propagation> function)

    caseTag = 'reflected';
    FIR_reflected = get_FIR_from_atm_TF( transferFunction_singleSideSpectrum_reflected, ray_delay_samples, fs, tag_auralization, caseTag );

end

%%% %% check plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P= 10.*log10( abs( transferFunction_singleSideSpectrum ).^2.);
% figure
% imagesc(P); colorbar('vert'); set(gca,'YDir','Normal')
% colormap('jet'); 
% caxis([-150 0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% apply atmospheric transfer function to noise emission using circular convolution

clear BlockOut;

% define synthesis block length [samples]
BlockLen = round(fs*dt_panam);

% convolute inputSignal with FIR of direct sound path
outputSignal_direct = overlapp_add_convolution( inputSignal, BlockLen, FIR_direct );

if considerGroundReflection == 1

    outputSignal_reflected = overlapp_add_convolution( inputSignal, BlockLen, FIR_reflected );
    outputSignal = outputSignal_direct + outputSignal_reflected;

else

    outputSignal = outputSignal_direct;

end

%%% %% check plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% plot(outputSignal);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if show == true
    tag_title =  ['OUTPUT - Spectrogram of auralized emission after propagation - ' tag_source];
    tag_save = ['_' tag_source '_spectrogram_emission_propagated'];
    PLOT_spectrogram(outputSignal, fs, tag_title, tag_auralization, tag_save);
else
end

% sound(outputSignal, 44100, 24 );

fprintf('\n- Applying propagation FIR filter/auralization time:\t%f sec\n', toc);

end

