function [L, t, P] = third_octave_level_history(x, fs, fc, T)
% THIRD_OCTAVE_LEVEL_HISTORY  Short-time 1/3-octave band levels (Leq over
% consecutive, non-overlapping frames of duration T). Used by the
% verification routines of the amplitude fluctuations.
%
%   [L, t, P] = third_octave_level_history(x, fs, fc, T)
%
% Band powers are obtained by summing the FFT power of each frame
% (rectangular window) over the bins between fc*2^(-1/6) and fc*2^(1/6).
% This is adequate for broadband test signals (e.g. white noise); it is not
% a standard-compliant filter bank.
%
% INPUTS
%   x  : [N x 1] signal
%   fs : sampling frequency [Hz]
%   fc : [1 x nB] band centre frequencies [Hz]
%   T  : frame duration [s], e.g. 0.1 (100 ms Leq, as in Lincke et al. 2025, Fig. 5)
%
% OUTPUTS
%   L  : [nFrames x nB] band levels [dB re. 1 (signal units)^2]
%   t  : [nFrames x 1] frame centre times [s]
%   P  : [nFrames x nB] band powers (linear)
%
% Author: Gil Felix Greco, Braunschweig (29.10.2026)
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to 
%   develop this code under the author's direction. All codes were verified 
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M  = round(T*fs);
nF = floor(numel(x)/M);
X  = reshape(x(1:nF*M), M, nF);
S  = abs(fft(X)).^2 / M^2;                   % [M x nF] two-sided power per bin
f  = (0:M-1).' * fs/M;

P = zeros(nF, numel(fc));
for b = 1:numel(fc)
    in = f >= fc(b)*2^(-1/6) & f < fc(b)*2^(1/6);
    P(:, b) = 2*sum(S(in, :), 1).';          % single-sided band power
end

L = 10*log10(P + realmin);
t = ((0:nF-1).' + 0.5)*T;

end
