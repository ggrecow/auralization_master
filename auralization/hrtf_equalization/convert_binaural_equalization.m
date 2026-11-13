function [y_out, h_conv] = convert_binaural_equalization( y_in, fs, eqFrom, eqTo )
% function [y_out, h_conv] = convert_binaural_equalization( y_in, fs, eqFrom, eqTo )
%
% Converts a binaural signal with one HRTF equalization (<eqFrom>) into the
% signal that would have been rendered with another one (<eqTo>), without
% re-rendering. Used by <MASTER_auralization_EngineAirframe> to obtain the
% diffuse-field ('DF') and free-field ('FF') equalized versions of the rendered
% (non-equalized, 'none') binaural signal (<generate_equalized_binaural> = 1 in
% the .ini file). It can also be applied to finished binaural .wav files, e.g.
% to compute metrics (LAmax, LAE, EPNL, ECMA-418-2 loudness, ...) from the FF
% version of the stimuli of a listening test:
%
%   y_ff = convert_binaural_equalization( y_bin, fs, 'none', 'FF' );
%
% Why this is exact (up to filter design tolerances):
%   - all equalizations ('none', 'DF', 'FF') are fixed, linear and time-invariant
%     filters per ear, applied in the same way to all time blocks and to the
%     direct and reflected paths;
%   - the block convolution of the framework (<overlapp_add_convolution>) is
%     linear and uses no output windowing, so a filter applied to every HRIR is
%     identical to the same filter applied to the output signal (checked with
%     <validate_hrtf_equalization>: < 0.005 dB per third-octave band).
%   The conversion filter per ear is the minimum-phase filter with magnitude
%   |E_to(f)| / |E_from(f)|, where E are the filters of <get_hrtf_equalization_filter>
%   (all of them minimum phase, so their phases follow from the magnitudes).
%   The binaural signals contain no energy above 22.05 kHz (the HRIRs are
%   resampled from 44.1 kHz), where the conversion gain is held constant.
%
% INPUTS
%   y_in   : [N x 2] binaural signal (left, right), in Pa
%   fs     : sampling frequency of y_in, Hz
%   eqFrom : equalization of y_in   ('none', 'DF' or 'FF')
%   eqTo   : desired equalization  ('none', 'DF' or 'FF')
%
% OUTPUTS
%   y_out  : [N x 2] converted binaural signal (same length as y_in)
%   h_conv : struct with the conversion filters at fs (.left, .right)
%
% Author: Gil Felix Greco, Braunschweig 11.11.2026
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to
%   develop this code under the author's direction. All codes were verified
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size( y_in, 2 ) ~= 2
    error( 'convert_binaural_equalization:input', 'y_in must be a [N x 2] binaural signal.' );
end

eqFrom = il_check_type( eqFrom );
eqTo   = il_check_type( eqTo );

if strcmp( eqFrom, eqTo )
    y_out  = y_in;
    h_conv = struct( 'left', 1, 'right', 1 );
    return
end

EQf = get_hrtf_equalization_filter( eqFrom );
EQt = get_hrtf_equalization_filter( eqTo );
fs_hrir = EQf.fs;

% design grid at the database fs, fine enough for the step-like FF filters
nTaps_hrir = max( [numel(EQf.left), numel(EQf.right), numel(EQt.left), numel(EQt.right), 256] );
Nd  = 2^nextpow2( 16*nTaps_hrir );
fd  = (0:Nd/2)' * fs_hrir/Nd;

% conversion filter length and design grid at fs
nTaps = 2^nextpow2( 2*nTaps_hrir*fs/fs_hrir );
N     = 2^nextpow2( 8*nTaps );
f     = (0:N/2)' * fs/N;

ears = {'left', 'right'};
y_out = zeros( size(y_in) );
for e = 1:2
    Af = abs( fft( EQf.(ears{e}), Nd ) );
    At = abs( fft( EQt.(ears{e}), Nd ) );
    R  = At(1:Nd/2+1) ./ max( Af(1:Nd/2+1), realmin );   % magnitude ratio at fs_hrir grid

    % to fs grid; constant above fs_hrir/2 (no signal energy there)
    Rf = interp1( fd, R, min( f, fd(end) ), 'linear' );

    h = il_min_phase( Rf, N, nTaps );
    h_conv.(ears{e}) = h;
    y_out(:,e) = fftfilt( h, y_in(:,e) );
end

end

%% inline function: minimum-phase FIR from a single-sided magnitude (real cepstrum)
function h = il_min_phase( mag, N, nTaps )
magFull = [mag; flipud( mag(2:end-1) )];
c       = real( ifft( log( max( magFull, realmin ) ) ) );
wc      = [1; 2*ones( N/2-1, 1 ); 1; zeros( N/2-1, 1 )];
h_min   = real( ifft( exp( fft( c .* wc ) ) ) );
h       = h_min( 1:nTaps );
nFade   = round( nTaps/10 );
fade    = 0.5 + 0.5*cos( pi*(1:nFade)'/nFade );
h(end-nFade+1:end) = h(end-nFade+1:end) .* fade;
end

%% inline function: check equalization type
function t = il_check_type( t )
switch upper( strtrim( t ) )
    case 'NONE', t = 'none';
    case 'DF',   t = 'DF';
    case 'FF',   t = 'FF';
    otherwise
        error( 'convert_binaural_equalization:type', ...
            'Unknown equalization ''%s'' (use ''none'', ''DF'' or ''FF'').', t );
end
end
