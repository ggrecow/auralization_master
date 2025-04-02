function outputSignal = apply_propagation(inputSignal, inputRayTracing, binaural_signal, show, tag_auralization, tag_source, considerGroundReflection )
% function outputSignal = apply_propagation(inputSignal, transferFunction, show, tag_auralization, tag_source, considerGroundReflection)
%
% The auralization of the aircraft noise at a receiver position is conducted here by applying the 
% atmospheric <transferFunction> computed using the fuction <get_propagation> to the <inputSignal> 
% (sinthesized sound source sound file). The transfer functions are
% transformed into FIR filter of a truncated length nTaps and applied to
% the <inputSignal> using uniformelly partitioned circular convolution.
%
% If required, generates stereo (binaural) signals. this feature uses measured HRTFs 
% from FABIAN database (zero-degree HATO) (DOI: 10.14279/depositonce-5718.5) 
%  to obtain
%  binaural signals from <inputSignal>. The aircraft trajectory is already indirectly accounted for  
% by the incidence angle directions, which were previously computed in the <get_propagation> 
% function, and are given as inputs here. Please refer to <AKhrirInterpolation>
% for the angle convention used by FABIAN (head-centered spherical
% coordinates). The HRTFs are applied individually to the direct and
% reflected ray paths.
%
% Head orientation is defined by changing the azimuth angles, by the
% <head_orientation> variable read in the input file.  By default, the front of the head 
%  points to positive x-axis, head_orientation = 0 (degrees). For ex.,
%   head_orientation = 90 -> head pointing to positive y-axis
%   head_orientation = 180 -> head pointing to negative x-axis
%   head_orientation = 270 -> head pointing to negative y-axis
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

global input_file

if binaural_signal == 1
    if isfield( input_file, 'head_orientation' )   
        head_orientation = str2double ( input_file.head_orientation );  % get <head_orientation> from <input_file>
    else
        head_orientation = 0; % default value
    end
end

tic;
%% Inputs

transferFunction = inputRayTracing.TF;
propagation_time = inputRayTracing.propagation_time;

% ray_delay = propagation_time(:,2) - propagation_time(:,1); % delay between reflected and direct ray paths, in seconds
ray_delay_samples = round ( ( propagation_time(:,2) .* fs) - ( propagation_time(:,1) .* fs) ); % delay between reflected and direct ray paths, in samples

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

% define synthesis block length [samples]
BlockLen = round(fs*dt_panam);

% convolute inputSignal with FIR of direct sound path
outputSignal_direct = overlapp_add_convolution( inputSignal, BlockLen, FIR_direct, 'mono' );

if binaural_signal == 1

    % incidence angle (spherical coordinates) of the incoming sound wave (direct path)
    azimuth = inputRayTracing.spherical_angles_HRTF.direct_path(:,1);
    elevation = inputRayTracing.spherical_angles_HRTF.direct_path(:,2);
    HATO = 0; % head-above-torso orientation, degrees

    % get HRIR
    [h.l_a, h.r_a] = AKhrirInterpolation(azimuth+head_orientation, elevation, HATO, 'measured_sh');

    % approach 1) overlap add
    outputSignal_direct_binaural =  overlapp_add_convolution( inputSignal, BlockLen, h, 'stereo' ); 

    clear h; h.l_a = FIR_direct; h.r_a = FIR_direct;
    outputSignal_direct_binaural =  overlapp_add_convolution( outputSignal_direct_binaural, BlockLen, h, 'stereo' );

    % approach 2) get binaural signal (direct sound path). Error: some artifacts, not sure why
    % load spherical harmonics data
    % d = load( fullfile( 'FABIAN_HRIR_measured_HATO_0' ) );
    % outputSignal_direct_binaural = AKshAura(outputSignal_direct, azimuth+head_orientation, elevation, d.SH.coeffLeft, d.SH.coeffRight, 'complex', true, 1, BlockLen);

end

if considerGroundReflection == 1

    % convolute inputSignal with FIR of reflected sound path
    outputSignal_reflected = overlapp_add_convolution( inputSignal, BlockLen, FIR_reflected, 'mono' );

    % output (mono) signal
    outputSignal.outputSignal = outputSignal_direct + outputSignal_reflected;

    if binaural_signal == 1

    % incidence angle (spherical coordinates) of the incoming sound wave (reflected path)
    clear azimuth elevation;
    azimuth = inputRayTracing.spherical_angles_HRTF.reflected_path(:,1);
    elevation = inputRayTracing.spherical_angles_HRTF.reflected_path(:,2);

    % get HRIRs
    clear h;
    [h.l_a, h.r_a] = AKhrirInterpolation(azimuth+head_orientation, elevation, HATO, 'measured_sh');

    % approach 1) overlap add
    outputSignal_reflected_binaural =  overlapp_add_convolution( inputSignal, BlockLen, h, 'stereo' ); 

    clear h; h.l_a = FIR_reflected; h.r_a = FIR_reflected;
    outputSignal_reflected_binaural =  overlapp_add_convolution( outputSignal_reflected_binaural, BlockLen, h, 'stereo' );

    % approach 2) get binaural signal (reflected sound path). Error: some artifacts, not sure why
    % outputSignal_reflected_binaural = AKshAura(outputSignal_reflected, azimuth+head_orientation, elevation, d.SH.coeffLeft, d.SH.coeffRight, 'complex', true, 1, BlockLen);
    
    % output (stereo) signal
    outputSignal.outputSignal_binaural = outputSignal_direct_binaural + outputSignal_reflected_binaural;

    end

else % only direct sound path is an output 

    outputSignal.outputSignal = outputSignal_direct;

    if binaural_signal == 1
        outputSignal.outputSignal_binaural = outputSignal_direct_binaural;
    end

end

%%% %% check plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% plot(outputSignal);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if show == true
    tag_title =  ['OUTPUT - Spectrogram of auralized emission after propagation - ' tag_source];
    tag_save = ['_' tag_source '_spectrogram_emission_propagated'];
    PLOT_spectrogram(outputSignal.outputSignal, fs, tag_title, tag_auralization, tag_save);
else
end

% sound(outputSignal.outputSignal, 44100, 24 );

if binaural_signal == 1
    if show == true
        tag_title =  ['OUTPUT - Spectrogram of auralized emission after propagation (binaural) - ' tag_source];
        tag_save = ['_' tag_source '_spectrogram_emission_propagated_binaural'];
        PLOT_spectrogram(outputSignal.outputSignal_binaural, fs, tag_title, tag_auralization, tag_save);
    else
    end
end

fprintf('\n- Applying propagation FIR filter/auralization time:\t%f sec\n', toc);

end

