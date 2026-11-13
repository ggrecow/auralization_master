function FF = get_hrtf_ff_equalization( opts )
% function FF = get_hrtf_ff_equalization( opts )
%
% Designs the free-field (FF) equalization filters of the FABIAN HRTFs: the
% inverted, minimum-phase FRONTAL HRTF (azimuth = 0 deg, elevation = 0 deg,
% head-centred coordinates of <AKhrirInterpolation>), one filter per ear.
%
% It is the free-field counterpart of the inverted CTF used for the
% diffuse-field (DF) equalization: same idea (invert, minimum phase), but
% with the frontal HRTF as reference instead of the diffuse-field average of
% all directions.
%
% Filtering the HRIR of direction (az, el) with these filters gives
%
%       H(az, el) / H(0, 0)          (per ear)
%
% i.e. the sound pressure of a FRONTAL plane wave, measured at the head
% position without the listener, that would produce the same ear signal.
% This is the input assumed by metrics that already contain the outer-ear
% effect of a frontal free field (e.g. ECMA-418-2 with field type
% 'free-frontal'), or that are defined on the free-field sound pressure
% without the listener (A-weighting, noy/PNL weighting). Same principle as
% ISO 11904-2 (conversion of ear signals to free-field related levels).
% For frontal incidence the ratio is 1, so the FF equalized binaural signal
% and the omnidirectional microphone give the same spectra and metrics.
%
% FF equalized signals are meant for COMPUTING METRICS. For headphone
% playback, use raw HRTFs + measured headphone filter, or DF equalized HRTFs
% with diffuse-field equalized headphones (see README, Sect. 1.4).
%
% This function is normally NOT called during the auralization: the filters
% are precomputed with utilities/generate_FABIAN_FF_equalization_mat.m and
% shipped in third_party/FABIAN_HRTF_DATABASE_v4/ (as for the CTF). It is
% called on the fly only if that .mat file is missing (see
% <get_hrtf_equalization_filter>). The FF equalized binaural signal is
% obtained from the rendered (non-equalized) signal with
% <convert_binaural_equalization>.
%
% INPUT (all fields optional)
%   opts.HATO         : head-above-torso orientation, deg            (default 0)
%   opts.fs_hrir      : sampling frequency of the FABIAN data, Hz    (default 44100)
%   opts.nfft         : FFT length used for the design               (default 16384)
%   opts.method       : how the frontal reference |H0| is inverted   (default 'band')
%       'band'  -> smooth filter whose gain is set per third-octave band, so that
%                  the band energy of the equalized frontal HRTF is exactly 0 dB
%                  in every band between f_lo and f_hi. For any other direction,
%                  the band levels equal those of |H(az,el)|^2/|H0|^2 averaged per
%                  band (band-wise FF referencing, as in ISO 11904-2). The pinna
%                  notches of H0 are NOT inverted into narrow peaks.
%       'exact' -> G = 1/|H0| on the FFT grid (no smoothing). Frontal incidence
%                  then gives a flat magnitude (exact BIN == MIC), but the deep
%                  notches of H0 become narrow peaks of G, which add energy for
%                  all other directions (for FABIAN, left ear: up to approx.
%                  +14 dB in the 6.3 kHz third-octave band). Use for checks only.
%       NOTE: a sliding 1/3-octave smoothing of |H0|^2 before inversion (as used
%       for the FABIAN CTF) was tested and rejected: for frontal incidence it
%       leaves band errors of up to 2.5 dB (8 kHz band), i.e. BIN ~= MIC.
%   opts.nIter        : max. iterations of the 'band' method          (default 50)
%   opts.f_lo         : lower limit of the equalization band, Hz     (default 50)
%   opts.f_hi         : upper limit of the equalization band, Hz     (default 16000)
%                       outside [f_lo, f_hi] the equalization fades
%                       to 0 dB over half an octave (raised cosine)
%   opts.max_boost_dB : maximum gain of the inverse filter, dB       (default 20)
%                       (limits the boost at deep pinna notches)
%   opts.nTaps        : length of the FIR filters, samples           (default 1024)
%
% OUTPUT
%   FF : struct with fields
%     .left, .right     : [nTaps x 1] minimum-phase FIR filters at fs_hrir
%     .fs               : sampling frequency of the filters (= opts.fs_hrir)
%     .f                : [nfft/2+1 x 1] frequency vector, Hz
%     .H0_left/.H0_right: smoothed magnitude of the frontal HRTF (linear)
%     .G_left/.G_right  : target magnitude of the filters (linear), i.e. the
%                         inverse of H0 with boost limit and band fades
%     .settings         : the options used (for documentation)
%
% REQUIREMENTS: AKtools (AKhrirInterpolation) and the FABIAN SH data
% (FABIAN_HRIR_measured_HATO_0.mat), both shipped with the framework.
%
% Author: Gil Felix Greco, Braunschweig 11.11.2026
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to
%   develop this code under the author's direction. All codes were verified
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% default options
if nargin < 1 || isempty( opts ), opts = struct(); end
defaults = struct( 'HATO', 0, 'fs_hrir', 44100, 'nfft', 16384, 'method', 'band', ...
                   'f_lo', 50, 'f_hi', 16000, 'max_boost_dB', 20, 'nTaps', 2048, 'nIter', 50, 'trans_oct', 0 );
names = fieldnames( defaults );
for k = 1:numel( names )
    if ~isfield( opts, names{k} ) || isempty( opts.(names{k}) )
        opts.(names{k}) = defaults.(names{k});
    end
end

if mod( opts.nfft, 2 ) ~= 0 || opts.nTaps > opts.nfft
    error( 'get_hrtf_ff_equalization:nfft', 'nfft must be even and >= nTaps.' );
end

%% frontal HRIRs, same spherical-harmonics representation as used in <get_FIR>
[h0L, h0R] = AKhrirInterpolation( 0, 0, opts.HATO, 'measured_sh' );

[FF.left,  FF.H0_left,  FF.G_left,  f] = il_design( h0L, opts );
[FF.right, FF.H0_right, FF.G_right]    = il_design( h0R, opts );

FF.fs       = opts.fs_hrir;
FF.f        = f;
FF.settings = opts;

end

%% inline function: reference magnitude, regularized inversion, minimum phase
function [h_eq, H0s, G, f] = il_design( h, opts )

N  = opts.nfft;
M  = N/2 + 1;                                  % single-sided bins
df = opts.fs_hrir / N;
f  = (0:M-1)' * df;

H = fft( h(:), N );
P = abs( H(1:M) ).^2;                          % power spectrum of the frontal HRTF

% fade weights: equalization fades to 0 dB outside [f_lo, f_hi] (half-octave raised cosine)
w  = ones( M, 1 );
t  = log2( max( f, realmin ) / opts.f_lo );   lo = t < 0;
w(lo) = 0.5 + 0.5*cos( pi * min( -t(lo), 0.5 ) / 0.5 );
t  = log2( max( f, realmin ) / opts.f_hi );   hi = t > 0;
w(hi) = 0.5 + 0.5*cos( pi * min(  t(hi), 0.5 ) / 0.5 );

switch lower( opts.method )

    case 'exact'
        % G = 1/|H0| on the FFT grid: exact narrowband cancellation of the frontal
        % HRTF. Deep pinna notches of H0 become narrow peaks of G, which add energy
        % for all other directions (up to approx. +14 dB in the 6.3 kHz third-octave
        % band for FABIAN). Use for the sanity check "frontal source: BIN == MIC".
        Ps   = P;
        G_dB = -10*log10( max( Ps, realmin ) );
        G_dB = min( G_dB, opts.max_boost_dB );

    case 'band'
        % G is a smooth curve through one gain per third-octave band, iterated so
        % that the BAND ENERGY of the equalized frontal HRTF is 0 dB in every band
        % between f_lo and f_hi. For other directions, the band levels of the
        % equalized HRTF then equal |H(az,el)|^2 / |H0|^2 averaged per band, i.e. the
        % band-wise free-field referencing of ISO 11904-2.
        fc   = 1000 * 2.^( (floor(3*log2(opts.f_lo/1000)) : ceil(3*log2(opts.f_hi/1000)))' / 3 );
        nB   = numel( fc );
        idxB = cell( nB, 1 );
        for b = 1:nB
            idxB{b} = find( f >= fc(b)*2^(-1/6) & f < fc(b)*2^(1/6) );
        end
        bandE = @(X) cellfun( @(id) mean( X(id) ), idxB );

        gB_dB = -10*log10( bandE( P ) );                 % first guess: inverse band energy
        for it = 1:opts.nIter
            G_dB  = il_interp_bands( fc, min( gB_dB, opts.max_boost_dB ), f, opts.trans_oct );
            res   = 10*log10( bandE( P .* 10.^( G_dB/10 ) ) );  % residual band level, dB
            gB_dB = gB_dB - res;
            if max( abs( res ) ) < 1e-3, break; end
        end
        G_dB = il_interp_bands( fc, min( gB_dB, opts.max_boost_dB ), f, opts.trans_oct );
        Ps   = 10.^( -G_dB/10 );                         % equivalent reference power

    otherwise
        error( 'get_hrtf_ff_equalization:method', ...
            'Unknown method ''%s'' (use ''band'' or ''exact'').', opts.method );
end
H0s = sqrt( Ps );

G  = 10.^( w .* G_dB / 20 );

% minimum-phase FIR (real cepstrum, same method as <get_diffuse_field_filter>)
magFull = [G; flipud( G(2:end-1) )];
c       = real( ifft( log( magFull ) ) );
wc      = [1; 2*ones( N/2-1, 1 ); 1; zeros( N/2-1, 1 )];
h_min   = real( ifft( exp( fft( c .* wc ) ) ) );

% truncate to nTaps, half-Hann fade-out over the last 10% of the taps
h_eq  = h_min( 1:opts.nTaps );
nFade = round( opts.nTaps/10 );
fade  = 0.5 + 0.5*cos( pi*(1:nFade)'/nFade );
h_eq(end-nFade+1:end) = h_eq(end-nFade+1:end) .* fade;

end

%% inline function: curve through third-octave band gains
function G_dB = il_interp_bands( fc, gB_dB, f, trans_oct )
% Piecewise-constant gain per third-octave band (band edges fc*2^(+-1/6)),
% with raised-cosine transitions of width <trans_oct> (octaves) centred on
% the band edges. Constant outside the first/last band (the band fades are
% applied afterwards).
x  = log2( max( f, realmin ) );
xe = log2( fc(1:end-1) ) + 1/6;              % log2 of the inner band edges
G_dB = gB_dB(1) * ones( size(x) );
for b = 1:numel( xe )
    step = gB_dB(b+1) - gB_dB(b);
    if trans_oct > 0
        u = min( max( (x - xe(b)) / trans_oct + 0.5, 0 ), 1 );   % 0 -> 1 across the edge
        s = 0.5 - 0.5*cos( pi*u );
    else
        s = double( x >= xe(b) );
    end
    G_dB = G_dB + step * s;
end
end
