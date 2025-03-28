function outputSignal = overlapp_add_convolution( inputSignal, BlockLen, IR )
% function outputSignal = overlapp_add_convolution( inputSignal, BlockLen, IR )
%
%   perform FFT-based block convolution using the overlap and add technique
%
%   inputs  
%   inputSignal : vector
%   
%   BlockLen : scalar
%   synthesis block length [samples]
%
%   IR : vector
%   impulse response, with size <nTaps>
%
%   output:
%   outputSignal : vector
%   convoluted output signal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get number of Blocks for convolution
L = ceil(size(inputSignal,1)/BlockLen);

% zero Pad input signal to integer divisor of block size N
inputSignal(end+1:BlockLen*L,:) = 0;

% get size of the FIR filter
nTaps = length(IR);

% allocate output signal
outputSignal = zeros(BlockLen*L+nTaps-1, 1);

% fft filt using overlap and add
for ll = 1:L

    outputSignal((ll-1)*BlockLen+1:ll*BlockLen+nTaps-1,:) = ...
                                                                                            outputSignal((ll-1)*BlockLen+1:ll*BlockLen+nTaps-1,:) ...
                                                                                            + fftfilt(inputSignal((ll-1)*BlockLen+1:ll*BlockLen,:), [IR(:,ll); zeros(BlockLen-1, 1)]);

end