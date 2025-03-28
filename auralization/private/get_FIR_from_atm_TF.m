function OUT = get_FIR_from_atm_TF( input, ray_delay_samples, fs, tag_auralization, caseTag )
% function FIR = get_FIR_from_atm_TF( input, ray_delay_samples, fs, tag_auralization, caseTag )
%
% get FIR filter from input, which is here an single-sided complex-valued
% atmospheric transfer function
%
%   1) input is single-sided complex-value spectrum --> get double-sided complex-valued spectrum
%   2) transform double-sided spectrum signal to time domain
%   3) shift impulse response to make it causal
%   4) truncate FIR to nTaps
%
%  nomenclature: 
%  nBins = number of frequency bins
%  nTimes : number of samples corresponding to time steps (or source/receiver combinations)
% 
%   inputs
%   input : vector [nBins x nTimes]
%           single-sided complex-value spectra
%
%   ray_delay_samples : vector
%          delay between reflected and direct ray paths, in samples
%
%   fs : scalar
%           sampling frequency
%
%   tag_auralization : string
%           name to save the plots
%
%   caseTag : string
%           flag to sinalize if we are dealing with 'direct' or 'reflected', so that we can save 
%           plots with correct name
%
%   output
%   FIR : vector
%       finite impulse response filter with size nTaps
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% lets define nTaps based on fs and a desired min. freq

min_freq = 20; % in (Hz)
nTaps = fs/20;

% check if nTaps is big enough
if max(ray_delay_samples) >=  nTaps
    nTapsNew = max(ray_delay_samples)*2; % arbitrary number
    warning('FIR filter: delay between reflected and direct rays is larger than nTaps = %.1g Adjusting nTaps to %.1g ', nTaps, nTapsNew);
    nTaps = nTapsNew; clear nTapsNew;
end

%% initialize inputs 

numTimeSteps = size(input,2);

% Preallocate output matrix for FIR filters
impulseResponse = zeros(size (input,1)*2-2 , numTimeSteps );

for i = 1:numTimeSteps

    % Get the transfer function for the current time step (column)
    H = input(:, i);

    % Compute full Hermitian symmetric spectrum
    H_full = [H; conj(flipud(H(2:end-1)))];

    % Compute IFFT to obtain impulse response
    impulseResponse(:, i)  = ifft(H_full, 'symmetric');

end

%% make it causal 
% https://www.youtube.com/watch?v=R6sKlHHTULM&list=LL&index=12&t=766s

% shift impulse response to make it simmetrical  
Nsamples = size(impulseResponse,1); % number of samples of the impulse response 

% sometimes it is necessary to mirror the IFFT(FRF), sometimes not. This loop
% tries to deal with that in a very coarse way. This problem/solution to be further
% investigated/improved
for i =1:numTimeSteps

    if impulseResponse(end,i) > 1e-8 % if this arbitrary threshold is achieved in the last bin, then it means we need to mirror the IR 
        impulseResponseShifted(:,i) = [impulseResponse(Nsamples/2+1:end,i); impulseResponse(1:Nsamples/2,i)  ];
    else
        impulseResponseShifted(:,i) = impulseResponse(:,i);
    end

end

   %% truncate FIR to nTaps

   % zero pad IR to avoid problems when windowing FIR to a desired number of taps
   impulseResponseShiftedPadded = [ zeros(nTaps, numTimeSteps); impulseResponseShifted; zeros(nTaps, numTimeSteps)];

   for i = 1:numTimeSteps

       [~, idxMax(i)] = max ( impulseResponseShiftedPadded(:,i) );
       impulseResponseCut(:,i) = impulseResponseShiftedPadded( idxMax(i) - (nTaps/2) : idxMax(i) + (nTaps/2) - 1, i); % ntaps/2

   end

   %% deal with phase differences
   
   switch caseTag
       case 'direct'

           for i = 1:numTimeSteps
               % minimum phase FIR
               OUT(:,i)  = AKphaseManipulation( impulseResponseCut(:,i) , fs, 'min', 0, 0);
           end

       case 'reflected'

           for i = 1:numTimeSteps
               % linear phase FIR
               OUT(:,i)  = AKphaseManipulation( impulseResponseCut(:,i) , fs, 'lin', ray_delay_samples(i), 0);
           end

   end

%% plot freq response of FIR filter

tBlock = round(numTimeSteps/2); % a time bin to plot the freq response of the FIR filter

switch caseTag
    case 'direct'

        % plot FIR fpr a fixed (single) source/receiver combination
        PLOT_FIR( OUT, input, tBlock, fs, [tag_auralization '_direct'] );

        % plot FIR spectrogram
        PLOT_FIR_spectrogram(OUT, input, fs, [tag_auralization '_direct'] );

    case 'reflected'

        % plot FIR fpr a fixed (single) source/receiver combination
        PLOT_FIR( OUT, input, tBlock, fs, [tag_auralization '_reflected'] );

        % plot FIR spectrogram
        PLOT_FIR_spectrogram(OUT, input, fs, [tag_auralization '_reflected'] );
end

end