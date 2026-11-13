% validate_hrtf_equalization.m
%
% Verification of the free-field (FF) equalization of the FABIAN HRTFs and of
% the conversion between equalizations (<convert_binaural_equalization>).
% The rendering chain of <get_FIR> + <overlapp_add_convolution> is reproduced
% here in a minimal form (HRIR -> equalization at 44.1 kHz -> resampling ->
% block convolution), without ray tracing.
%
%   1) frontal HRTF x FF filter: third-octave band levels = 0 dB (per ear)
%   2) FF equalized HRTFs for several directions vs. band-wise referencing
%      |H(az,el)|^2/|H0|^2 per third-octave band (ISO 11904-2 principle)
%   3) frontal source: FF equalized binaural signal vs. mono (microphone):
%      third-octave band levels and LAeq must match
%   4) time-varying flyover chain: FF applied to the HRIRs while rendering
%      (reference) vs. FF obtained afterwards from the 'none' and 'DF'
%      renderings with <convert_binaural_equalization>, as done by the framework
%   5) figure: frontal HRTF, DF and FF equalization filters (optional)
%
% Results obtained when the filters were generated (Octave, 11.11.2026), all
% third-octave bands 50 Hz - 12.5 kHz: 1) max |dL| = 0.02 dB; 2) max |dL| = 0.11 dB;
% 3) max |dL| = 0.08 dB per band, |dLAeq| < 0.07 dB; 4) max |dL| = 0.002 dB per
% band, |dLAeq| < 0.001 dB. Residuals come from the finite FFT resolution at the
% band edges and the truncation of the FF filters to 2048 taps.
%
% Author: Gil Felix Greco, Braunschweig 11.11.2026
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to
%   develop this code under the author's direction. All codes were verified
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

do_plot = 1;        % 1: plot the filters (step 5)
fs      = 48000;    % framework sampling frequency used in steps 3 and 4
fs_hrir = 44100;    % FABIAN

EQ.none = get_hrtf_equalization_filter( 'none' );
EQ.DF   = get_hrtf_equalization_filter( 'DF' );
EQ.FF   = get_hrtf_equalization_filter( 'FF' );
fprintf( 'FF filters: %s (%d taps)\n', EQ.FF.source, numel(EQ.FF.left) );

fc  = 1000 * 2.^( (-13:11)/3 );     % third-octave bands 50 Hz - 12.5 kHz
lev = @(P, f) arrayfun( @(k) 10*log10( sum( P( f >= fc(k)*2^(-1/6) & f < fc(k)*2^(1/6) ) ) ), 1:numel(fc) );

%% 1) frontal HRTF x FF filter
N  = 2^16;
f  = (0:N/2)' * fs_hrir/N;
[h0L, h0R] = AKhrirInterpolation( 0, 0, 0, 'measured_sh' );
P0L = abs( fft( h0L, N ) ).^2;           P0L = P0L(1:N/2+1);
P0R = abs( fft( h0R, N ) ).^2;           P0R = P0R(1:N/2+1);
EL  = abs( fft( EQ.FF.left,  N ) ).^2;   EL  = EL(1:N/2+1);
ER  = abs( fft( EQ.FF.right, N ) ).^2;   ER  = ER(1:N/2+1);
d1  = [ lev( P0L.*EL, f ) - lev( ones(size(f)), f ); lev( P0R.*ER, f ) - lev( ones(size(f)), f ) ];
fprintf( '\n1) frontal HRTF x FF filter, band levels (L, R): max |dL| = %.3f dB\n', max( abs( d1(:) ) ) );

%% 2) other directions: FF equalized HRTF vs. band-wise referencing
az = [ 0  0 90 180  0 45 270];
el = [90 45  0   0 -30 30  0];
[hL, hR] = AKhrirInterpolation( az, el, 0, 'measured_sh' );
PL = abs( fft( hL, N ) ).^2;  PL = PL(1:N/2+1,:);
PR = abs( fft( hR, N ) ).^2;  PR = PR(1:N/2+1,:);
fprintf( '\n2) FF equalized HRTF vs. band-wise referencing |H|^2/|H0|^2 per band\n' );
for k = 1:numel(az)
    flat = lev( ones(size(f)), f );
    dL = lev( PL(:,k).*EL, f ) - flat - ( lev( PL(:,k), f ) - lev( P0L, f ) );
    dR = lev( PR(:,k).*ER, f ) - flat - ( lev( PR(:,k), f ) - lev( P0R, f ) );
    fprintf( '   az = %4d, el = %4d: max |dL| = %.3f dB\n', az(k), el(k), max( abs( [dL dR] ) ) );
end

%% 3) and 4): minimal rendering chain
dt = 0.1;  nT = 30;  BL = round( fs*dt );
rng( 1 );  x = randn( BL*nT, 1 );            % white noise source signal

% 3) frontal source, direct path only, constant delay
hrir = @(eq, a, e) il_hrir( a, e, EQ.(eq), fs, fs_hrir );
[bL, bR] = hrir( 'FF', 0, 0 );
y_ff = [ conv( x, bL ), conv( x, bR ) ];
y_mono = [ x; zeros( size(y_ff,1) - numel(x), 1 ) ];
fy = (0:size(y_ff,1)-1)' * fs/size(y_ff,1);
Yff = abs( fft( y_ff ) ).^2;  Ym = abs( fft( y_mono ) ).^2;
d3  = [ lev( Yff(:,1), fy ) - lev( Ym, fy ); lev( Yff(:,2), fy ) - lev( Ym, fy ) ];
dLA = [ il_LA( y_ff(:,1), fs ) il_LA( y_ff(:,2), fs ) ] - il_LA( y_mono, fs );
fprintf( '\n3) frontal source, FF binaural vs. mono: max |dL| per band = %.3f dB, dLAeq (L, R) = %+.3f, %+.3f dB\n', ...
    max( abs( d3(:) ) ), dLA );

% 4) flyover: front -> overhead -> back, direct path, time-varying HRIRs
t   = linspace( -1, 1, nT );
azT = 180*( t > 0 );   elT = 90 - 80*abs( t );
types = {'none', 'DF', 'FF'};
for m = 1:3
    y.(types{m}) = zeros( BL*nT + 4096*4, 2 );
    for k = 1:nT
        [bL, bR] = hrir( types{m}, azT(k), elT(k) );
        idx = (k-1)*BL + ( 1:BL + numel(bL) - 1 );
        xk  = x( (k-1)*BL + (1:BL) );
        y.(types{m})(idx,:) = y.(types{m})(idx,:) + [ conv( xk, bL ), conv( xk, bR ) ];
    end
end
y.none2FF = convert_binaural_equalization( y.none, fs, 'none', 'FF' );
y.DF2FF   = convert_binaural_equalization( y.DF,   fs, 'DF',   'FF' );
fprintf( '\n4) flyover, FF while rendering vs. FF by post-filter conversion\n' );
for src = {'none2FF', 'DF2FF'}
    a = y.(src{1});  b = y.FF;
    fa = (0:size(a,1)-1)' * fs/size(a,1);
    A = abs( fft( a ) ).^2;  B = abs( fft( b ) ).^2;
    d4 = [ lev( A(:,1), fa ) - lev( B(:,1), fa ); lev( A(:,2), fa ) - lev( B(:,2), fa ) ];
    dLA4 = [ il_LA( a(:,1), fs ) - il_LA( b(:,1), fs ), il_LA( a(:,2), fs ) - il_LA( b(:,2), fs ) ];
    fprintf( '   %-7s: max |dL| per band = %.3f dB, dLAeq (L, R) = %+.4f, %+.4f dB\n', src{1}, max( abs( d4(:) ) ), dLA4 );
end
fprintf( '   LAeq binaural - LAeq mono (L, R): none %+.1f %+.1f | DF %+.1f %+.1f | FF %+.1f %+.1f dB (white noise)\n', ...
    [ il_LA( y.none, fs ); il_LA( y.DF, fs ); il_LA( y.FF, fs ) ]' - il_LA( [x; zeros(1e4,1)], fs ) );

%% 5) figure
if do_plot
    D  = abs( fft( EQ.DF.left, N ) );  D = D(1:N/2+1);
    figure( 'Name', 'FABIAN HRTF equalization' );
    semilogx( f, 10*log10(P0L), f, 10*log10(P0R), f, 20*log10(D), '--', ...
              f, 10*log10(EL), ':', f, 10*log10(ER), ':' );
    grid on; xlim( [50 20000] ); ylim( [-25 25] );
    xlabel( 'Frequency (Hz)' ); ylabel( 'Magnitude (dB)' );
    legend( 'frontal HRTF, left', 'frontal HRTF, right', 'DF filter (inverted CTF)', ...
            'FF filter, left', 'FF filter, right', 'Location', 'southwest' );
    title( 'FABIAN (HATO 0): frontal HRTF and equalization filters' );
end

%% local functions
function [bL, bR] = il_hrir( az, el, EQ, fs, fs_hrir )
% HRIR processing as in <get_FIR> (interpolation, resampling with preserved
% frequency response), plus an equalization of the HRIRs (full linear
% convolution) to build the reference renderings
[bL, bR] = AKhrirInterpolation( az, el, 0, 'measured_sh' );
bL = conv( bL, EQ.left );  bR = conv( bR, EQ.right );
if fs ~= fs_hrir
    bL = resample( bL, fs, fs_hrir ) * fs_hrir/fs;
    bR = resample( bR, fs, fs_hrir ) * fs_hrir/fs;
end
end

function LA = il_LA( y, fs )
% A-weighted level (dB re 1, energy over the whole signal), frequency domain
N  = size( y, 1 );
f  = (0:N-1)' * fs/N;  f = min( f, fs - f );
RA = @(f) ( 12194^2 * f.^4 ) ./ ( (f.^2 + 20.6^2) .* sqrt( (f.^2 + 107.7^2) .* (f.^2 + 737.9^2) ) .* (f.^2 + 12194^2) );
W  = ( RA(f) / RA(1000) ).^2;
LA = 10*log10( sum( abs( fft( y ) ).^2 .* W, 1 ) );
end
