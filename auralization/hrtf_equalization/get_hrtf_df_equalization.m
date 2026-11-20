function DF = get_hrtf_df_equalization( hrtf, opts )
% function DF = get_hrtf_df_equalization( hrtf, opts )
%
% Designs the diffuse-field (DF) equalization filter of an HRTF set: the inverted,
% 3rd-octave smoothed, minimum-phase common transfer function (CTF), i.e. the
% diffuse-field average of the HRTFs over all directions and both ears. Same filter for
% both ears, as the CTF provided with the FABIAN database.
%
% The diffuse-field power is obtained directly from the spherical-harmonics (SH)
% coefficients of the measured HRTFs. For orthonormal SH basis functions (AKtools,
% 'complex'), the average of |H|^2 over the sphere is (Parseval)
%
%       P_ear(f) = 1/(4*pi) * sum_{n,m} |f_nm(f)|^2,          P_DF = (P_left + P_right)/2
%
% so no sampling grid or quadrature weights are needed. The magnitude of the HRTFs is not
% changed by the time alignment of the HUTUBS data before the SH transform.
% Steps: P_DF on the 129 frequencies of the SH data -> linear interpolation to the design
% grid -> sliding third-octave smoothing (power average over f*2^(+-1/6)) -> inversion
% (gain limited to <max_boost_dB>) -> minimum phase (real cepstrum) -> <nTaps> samples.
%
% Used for hrtf_database = HUTUBS, for which no CTF is provided (see
% <get_hrtf_equalization_filter>). For FABIAN, the framework uses the CTF provided with the
% database; this function can also be called for FABIAN, to check that the method
% reproduces that CTF (see <validate_hrtf_database>: within 0.2 dB per third-octave band
% between 50 Hz and 16 kHz).
%
% INPUTS
%   hrtf : struct from <get_hrtf_database>
%   opts : struct (all fields optional)
%     .nfft         : FFT length used for the design          (default 4096)
%     .nTaps        : length of the FIR filter, samples        (default 512)
%     .max_boost_dB : maximum gain of the filter, dB           (default 20)
%
% OUTPUT
%   DF : struct with fields
%     .left, .right : [nTaps x 1] minimum-phase FIR filter at DF.fs (identical)
%     .fs           : sampling frequency of the filter (= hrtf.fs)
%     .f            : [nfft/2+1 x 1] frequency vector of the design grid, Hz
%     .P_df         : [nfft/2+1 x 1] smoothed diffuse-field power (linear, re 1)
%     .G            : [nfft/2+1 x 1] target magnitude of the filter (linear)
%     .settings     : options used, and the HRTF set (.hrtf = hrtf.name)
%
% Gil Felix Greco, Braunschweig 19.11.2026
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to develop this
%   code under the author's direction. All codes were verified by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2 || isempty( opts ), opts = struct(); end
defaults = struct( 'nfft', 4096, 'nTaps', 512, 'max_boost_dB', 20 );
names = fieldnames( defaults );
for k = 1:numel( names )
    if ~isfield( opts, names{k} ) || isempty( opts.(names{k}) )
        opts.(names{k}) = defaults.(names{k});
    end
end
if mod( opts.nfft, 2 ) ~= 0 || opts.nTaps > opts.nfft
    error( 'get_hrtf_df_equalization:nfft', 'nfft must be even and >= nTaps.' );
end

%% SH coefficients of the measured HRTFs, [(N+1)^2 x nFreq] per ear

switch hrtf.database
    case 'HUTUBS'
        cL  = hrtf.sh.HRIR.CoeffLeft;
        cR  = hrtf.sh.HRIR.CoeffRight;
        fSH = double( hrtf.sh.HRIR.f(:) );
    case 'FABIAN'
        S   = load( hrtf.file, 'SH' );
        cL  = S.SH.coeffLeft;
        cR  = S.SH.coeffRight;
        fSH = double( S.SH.f(:) );
    otherwise
        error( 'get_hrtf_df_equalization:database', 'Unknown HRTF database ''%s''.', hrtf.database );
end

%% diffuse-field power (Parseval, orthonormal SH)

P_sh = ( sum( abs( cL ).^2, 1 ) + sum( abs( cR ).^2, 1 ) ).' / ( 2*4*pi );   % [nFreq x 1]

%% design grid, third-octave smoothing, inversion

N  = opts.nfft;
f  = (0:N/2)' * hrtf.fs/N;
Pf = interp1( fSH, P_sh, f, 'linear', P_sh(end) );     % constant above the last SH frequency

P_df = Pf;
for k = 2:numel( f )                                   % 0 Hz bin kept
    id = f >= f(k)*2^(-1/6) & f <= f(k)*2^(1/6);
    P_df(k) = mean( Pf(id) );
end

G = min( 1 ./ sqrt( max( P_df, realmin ) ), 10^( opts.max_boost_dB/20 ) );

%% minimum-phase FIR (real cepstrum), truncated with a raised-cosine fade-out

magFull = [G; flipud( G(2:end-1) )];
c       = real( ifft( log( magFull ) ) );
wc      = [1; 2*ones( N/2-1, 1 ); 1; zeros( N/2-1, 1 )];
h_min   = real( ifft( exp( fft( c .* wc ) ) ) );

h     = h_min( 1:opts.nTaps );
nFade = round( opts.nTaps/10 );
fade  = 0.5 + 0.5*cos( pi*(1:nFade)'/nFade );
h(end-nFade+1:end) = h(end-nFade+1:end) .* fade;

%% output

DF.left          = h;
DF.right         = h;
DF.fs            = hrtf.fs;
DF.f             = f;
DF.P_df          = P_df;
DF.G             = G;
DF.settings      = opts;
DF.settings.hrtf = hrtf.name;

end
