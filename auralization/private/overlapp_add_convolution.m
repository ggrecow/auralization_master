function outputSignal = overlapp_add_convolution( inputSignal, BlockLen, IR, type, irIndex )
% function outputSignal = overlapp_add_convolution( inputSignal, BlockLen, IR, type, irIndex )
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
%   type  string
%   'mono' or 'stereo'
%
%   irIndex : vector (OPTIONAL)
%   column of <IR> used for each block, i.e. block ll is filtered with IR(:,irIndex(ll)).
%   DEFAULT: irIndex = 1:nBlocks (block ll uses column ll). Used with
%   <doppler_method> = delay_line, where the blocks are in reception time and each
%   block uses the filter of the emission time step heard in it (see <apply_propagation>)
%
%   output:
%   outputSignal : vector
%   convoluted output signal
%
% Gil Felix Greco, Braunschweig 16.11.2026 - optional input <irIndex>
%
%   AI disclosure (16.11.2026 changes): Generative AI assistance (Anthropic Claude) was
%   used to develop this code under the author's direction. All codes were verified by
%   the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get number of Blocks for convolution
L = floor(size(inputSignal,1)/BlockLen);

% filter (column of IR) used for each block
if nargin < 5 || isempty(irIndex)
    irIndex = 1:L;
end
if numel(irIndex) < L
    error( 'auralization:irIndexSize', ...
        'overlapp_add_convolution: <irIndex> has %d elements, but %d blocks are processed.', numel(irIndex), L );
end

% zero Pad input signal to integer divisor of block size N
inputSignal(end+1:BlockLen*L,:) = 0;

switch type
    case 'mono'

        % get size of the FIR filter
        nTaps = length(IR);

        % allocate output signal
        outputSignal = zeros(BlockLen*L+nTaps-1, 1);

        % fft filt using overlap and add
        for ll = 1:L

            outputSignal((ll-1)*BlockLen+1:ll*BlockLen+nTaps-1,:) = outputSignal((ll-1)*BlockLen+1:ll*BlockLen+nTaps-1,:) ...
                                                                    + fftfilt(inputSignal((ll-1)*BlockLen+1:ll*BlockLen,:), [IR(:,irIndex(ll)); zeros(BlockLen-1, 1)]);

        end

    case 'stereo'

        % get size of the FIR filter
        nTaps = size(IR.leftEar,1);

        % allocate output signal
        outputSignal = zeros(BlockLen*L+nTaps-1, 2);

        % fft filt using overlap and add
        for ll = 1:L

            outputSignal((ll-1)*BlockLen+1:ll*BlockLen+nTaps-1,:) = outputSignal((ll-1)*BlockLen+1:ll*BlockLen+nTaps-1,:) ...
                                                                    + fftfilt(inputSignal((ll-1)*BlockLen+1:ll*BlockLen,:), [[IR.leftEar(:,irIndex(ll)) IR.rightEar(:,irIndex(ll))]; zeros(BlockLen-1, 2)]);

        end

end % end <overlapp_add_convolution> function





