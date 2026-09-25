function [h_df, H_df, f] = get_diffuse_field_filter( fs, nGrid, fracOct, f_lim, nTaps )
% function [h_df, H_df, f] = get_diffuse_field_filter( fs, nGrid, fracOct, f_lim, nTaps )
%
% Computes a minimum-phase diffuse-field (DF) equalization filter for the
% FABIAN HRIRs (HATO = 0), to be applied to ALL HRIRs (both ears, all
% directions). The resulting directional transfer functions (DTF = HRTF/DF)
% keep all direction-dependent cues (ITD, ILD, pinna notches) unchanged,
% since the same filter is applied everywhere; only the direction-
% independent part (e.g. ear-canal / concha resonance) is removed.
%
% Procedure:
%   1) HRIRs on a quasi-uniform (equal-area Fibonacci) grid over the sphere
%   2) resample to the framework sampling frequency <fs> (FABIAN: 44.1 kHz)
%   3) DF magnitude = RMS over all directions and both ears
%   4) fractional-octave smoothing of the DF power spectrum
%   5) inversion, with the inverse held constant outside f_lim (no boost of
%      the low-frequency region or of the high-frequency roll-off)
%   6) minimum-phase inverse filter (real cepstrum), truncated to nTaps
%
% INPUTS (all optional):
%   fs      : sampling frequency of the framework [Hz]          (default 48000)
%   nGrid   : number of directions used for the DF average      (default 2000)
%   fracOct : smoothing bandwidth, 1/fracOct octave             (default 3)
%   f_lim   : [f_low f_high], band where the inverse is applied (default [200 16000]) Hz
%   nTaps   : length of the output FIR filter                   (default 512)
%
% OUTPUTS:
%   h_df : [nTaps x 1] minimum-phase DF equalization filter (at fs)
%   H_df : smoothed DF magnitude (before inversion), single-sided
%   f    : frequency vector of H_df [Hz]
%
% The filter is cached (persistent) for identical inputs, since it only
% depends on the HRTF database and not on the auralization case.
%
% Gil Felix Greco, Braunschweig 25.09.2026 
%
% AI disclosure: code development assisted 
% by Claude Opus 5 (Anthropic). All codes were verified by 
% the authors.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1 || isempty(fs),      fs      = 48000;       end
if nargin < 2 || isempty(nGrid),   nGrid   = 2000;        end
if nargin < 3 || isempty(fracOct), fracOct = 3;           end
if nargin < 4 || isempty(f_lim),   f_lim   = [200 16000]; end
if nargin < 5 || isempty(nTaps),   nTaps   = 512;         end

persistent cache
key = sprintf('%g_%d_%g_%g_%g_%d', fs, nGrid, fracOct, f_lim(1), f_lim(2), nTaps);
if isstruct(cache) && isfield(cache, 'key') && strcmp(cache.key, key)
    h_df = cache.h_df; H_df = cache.H_df; f = cache.f;
    return
end

fs_hrir = 44100; % sampling frequency of the FABIAN database

%% 1) equal-area Fibonacci grid (uniform weights -> plain mean is a surface average)
k  = (0:nGrid-1)';
el = asind( 1 - 2*(k+0.5)/nGrid );     % elevation [deg], 90 = north pole, 0 = front
az = mod( k*137.50776405, 360 );       % azimuth [deg], golden-angle increments

[l, r] = AKhrirInterpolation( az, el, 0, 'measured_sh' );

%% 2) resample to framework fs (same processing as the HRIRs used in get_FIR)
if fs ~= fs_hrir
    l = resample( l, fs, fs_hrir );
    r = resample( r, fs, fs_hrir );
end

%% 3) diffuse-field power average (all directions, both ears)
nFFT = 2^nextpow2( max( 4*size(l,1), 2*nTaps ) );
Hl   = fft( l, nFFT );
Hr   = fft( r, nFFT );
P    = mean( [abs(Hl).^2, abs(Hr).^2], 2 );

nSingle = nFFT/2 + 1;
P = P(1:nSingle);
f = (0:nSingle-1)' * fs / nFFT;

%% 4) fractional-octave smoothing (on power)
Ps = P;
bw = 2^( 1/(2*fracOct) );
for i = 2:nSingle
    idx   = f >= f(i)/bw & f <= f(i)*bw;
    Ps(i) = mean( P(idx) );
end
H_df = sqrt( Ps );

%% 5) inversion, held constant outside f_lim
Hinv = 1 ./ H_df;
iLow  = find( f >= f_lim(1), 1, 'first' );
iHigh = find( f <= f_lim(2), 1, 'last'  );
Hinv(1:iLow-1)    = Hinv(iLow);
Hinv(iHigh+1:end) = Hinv(iHigh);

%% 6) minimum-phase FIR via real cepstrum
magFull = [Hinv; flipud( Hinv(2:end-1) )];            % double-sided magnitude
c       = ifft( log( max(magFull, eps) ), 'symmetric' );  % real cepstrum
w       = [1; 2*ones(nFFT/2-1, 1); 1; zeros(nFFT/2-1, 1)]; % fold to causal part
h_min   = real( ifft( exp( fft( c .* w ) ) ) );

% truncate with half-Hann fade-out over the last 25% of the taps
h_df    = h_min(1:nTaps);
nFade   = round( nTaps/4 );
fade    = hann( 2*nFade );
h_df(end-nFade+1:end) = h_df(end-nFade+1:end) .* fade(nFade+1:end);

cache = struct( 'key', key, 'h_df', h_df, 'H_df', H_df, 'f', f );

end