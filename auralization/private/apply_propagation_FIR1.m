function outputSignal = apply_propagation_FIR1(inputSignal, transferFunction, nTaps, show, tag_auralization, tag_source, considerGroundReflection )
% function outputSignal = apply_propagation_FIR1(inputSignal, transferFunction, nTaps, show, tag_auralization, tag_source, considerGroundReflection)
%
% The auralization of the aircraft noise at a receiver position is conducted here by applying the 
% atmospheric <transferFunction> computed using the fuction <get_propagation> to the <inputSignal> 
% (sinthesized sound source sound file). The transfer functions are
% transformed into FIR filter of a truncated length nTaps and applied to
% the <inputSignal> using uniformelly partitioned circular convolution.
%
%
% INPUT:
%       inputSignal : [N x 1] vector
%       Auralized source description, in Pascal 
%
%       transferFunction : cell
%       results from ART (ray-tracing) propagation tool containing the
%       atmospheric transfer fuction for each Source x receiver
%       combinations. For each combination, the transfer fuction of the direct ray (1st column), 
%       1st order reflection ray (2nd column) and their combination (3rd column) is provided  
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global fs % sampling frequency used for auralization, defined in the <auralization_master>
global dt_panam % dt from panam, defined in the <auralization_master>

tic;
%% Get atmospheric transfer function 

% declare variables to pre allocate memory
nTF = size(transferFunction,1);
TF_direct = zeros( size(transferFunction,1), size(transferFunction{1,1}.freq,1) );
% TF_reflection = zeros( size(transferFunction,1), size(transferFunction{1,1}.freq,1) );
TF_combined = zeros( size(transferFunction,1), size(transferFunction{1,1}.freq,1) );

for i = 1:nTF
    TF_direct(i,:) = transferFunction{i,1}.freq(:,1);
    % TF_reflection(i,:) = transferFunction{i,1}.freq(:,2);
    TF_combined(i,:) = transferFunction{i,1}.freq(:,3);
end

% freq vector from FR - single sided spectrum - nBins = (round(fs*dt_panam)+1)/2 (see <get_propagation> function)
freq = transferFunction{1,1}.freqVector(:,1);

% number of blocks (i.e. time instances)
nTime = size(transferFunction,1);  

clear nTF

%% prepare atmospheric transfer function for circular convolution

% get atmospheric transfer function (single side spectrum)

if considerGroundReflection == 0

% direct sound 
transferFunction_singleSideSpectrum = TF_direct; % single sided spectrum - nBins = (round(fs*dt_panam)+1)/2 (see <get_propagation> function)

elseif considerGroundReflection == 1

% direct sound + 1st order reflection  
transferFunction_singleSideSpectrum = TF_combined; % single sided spectrum - nBins = (round(fs*dt_panam)+1)/2 (see <get_propagation> function)

end

% lets work only with row vectors
transferFunction_singleSideSpectrum = transpose(transferFunction_singleSideSpectrum);

%%% %% check plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P= 10.*log10( abs( transferFunction_singleSideSpectrum ).^2.);
% figure
% imagesc(P); colorbar('vert'); set(gca,'YDir','Normal')
% colormap('jet'); 
% caxis([-150 0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% double side the atmospheric transfer function  
for i = 1:size(transferFunction_singleSideSpectrum,2) % loop over time steps
    transferFunction_doubleSideSpectrum(:,i) =  [transferFunction_singleSideSpectrum(:,i); conj(flipud(transferFunction_singleSideSpectrum(2:end-1,i)))];
end

% transferFunction_doubleSideSpectrum = real(transferFunction_doubleSideSpectrum);

%%% %% check plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P= 10.*log10( abs( transferFunction_doubleSideSpectrum ).^2.);
% figure
% imagesc(P); colorbar('vert'); set(gca,'YDir','Normal')
% colormap('jet'); 
% caxis([-150 0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% transform double side spectrum signal to time domain
impulseResponse = ifft(transferFunction_doubleSideSpectrum, 'symmetric');

% shift impulse response to make it simmetrical  
Nsamples = size(impulseResponse,1); % number of samples of the impulse response 
% impulseResponseShifted = [impulseResponse(Nsamples/2+1:end,:); impulseResponse(1:Nsamples/2,:)  ];

% sometimes it is necessary to mirror the IFFT(FRF), sometimes not. This loop
% tries to deal with that in a very coarse way. This problem/solution to be further
% investigated/improved
for i =1:nTime

    if impulseResponse(end,i) > 1e-8 % if this arbitrary threshold is achieved in the last bin, then it means we need to mirror the IR 
        impulseResponseShifted(:,i) = [impulseResponse(Nsamples/2+1:end,i); impulseResponse(1:Nsamples/2,i)  ];
    else
        impulseResponseShifted(:,i) = impulseResponse(:,i);
    end

end

%% normalize responses directly from transfer function

clear impulseResponseCut

% nTaps = 2^14; 
% nTaps = round(fs*dt_panam);  

% zero pad shifted IR to avoid problems when windowing FIR to a desired number of taps 
impulseResponseShiftedPadded = [ zeros(nTaps, nTime); impulseResponseShifted; zeros(nTaps, nTime)];

for i = 1:size(impulseResponseShiftedPadded,2)

    [~, idxMax(i)] = max ( impulseResponseShiftedPadded(:,i) );
    impulseResponseCut(:,i) = impulseResponseShiftedPadded( idxMax(i) - (nTaps/2) : idxMax(i) + (nTaps/2) - 1, i); % ntaps/2

end

% window IR to get zeros on the edges
% impulseResponseCut = impulseResponseCut.*hann(size(impulseResponseCut,1),'symmetric');  

%%% %% check plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure( 'name', 'original impulseResponse'), plot(impulseResponse(:,21)); % plot original impulseResponse
% figure( 'name', 'impulseResponseShifted'), plot(impulseResponseShifted(:,21)); % plot shifted impulseResponse
% figure( 'name', 'impulseResponseShiftedPadded'),  plot(impulseResponseShiftedPadded); % plot shifted impulseResponse
% figure( 'name', 'impulseResponseCut'), plot(impulseResponseCut(:,21)); % plot impulseResponseCut
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot freq response of FIR filter

% plot FIR fpr a fixed (single) source/receiver combination
tBlock = round(length(transferFunction)/2); % a time bin to plot the freq response of the FIR filter
PLOT_FIR( impulseResponseCut, transferFunction_singleSideSpectrum, tBlock, fs, tag_auralization );

% plot FIR spectrogram
PLOT_FIR_spectrogram(impulseResponseCut, transferFunction_singleSideSpectrum, fs, tag_auralization);

%% apply atmospheric transfer function to noise emission using circular convolution

clear BlockOut;

% define synthesis block length [samples]
BlockLen = round(fs*dt_panam);               

% get number of Blocks for convolution
L = ceil(size(inputSignal,1)/BlockLen);

% zero Pad input signal to integer divisor of block size N
inputSignal(end+1:BlockLen*L,:) = 0;

% allocate output signal
outputSignal = zeros(BlockLen*L+nTaps-1, 1);

% Impulse response
IR = impulseResponseCut;

% fft filt using overlap and add
for ll = 1:L

    outputSignal((ll-1)*BlockLen+1:ll*BlockLen+nTaps-1,:) = ...
                                                                                            outputSignal((ll-1)*BlockLen+1:ll*BlockLen+nTaps-1,:) ...
                                                                                            + fftfilt(inputSignal((ll-1)*BlockLen+1:ll*BlockLen,:), [IR(:,ll); zeros(BlockLen-1, 1)]);

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

