% generate_FABIAN_FF_equalization_mat.m
%
% MAINTAINER SCRIPT (not needed by users): designs the free-field (FF)
% equalization filters of the FABIAN HRTFs (inverted frontal HRTF, one filter
% per ear, see <get_hrtf_ff_equalization>) and saves them to
%
%   third_party/FABIAN_HRTF_DATABASE_v4/FABIAN_FF_equalization_measured_HATO_0.mat
%
% The .mat file is shipped with the framework (as the FABIAN CTF used for the
% diffuse-field equalization), so that:
%   - the filters used in all auralizations are fixed and documented
%     (reproducibility; results do not depend on design code changes);
%   - they can be used in post-processing scripts without AKtools/FABIAN
%     (e.g. <convert_binaural_equalization>, metric scripts, figures).
% Designing the filters on the fly takes < 1 s, so speed is not the reason.
%
% Run this script again (and commit the .mat) whenever the design settings in
% <get_hrtf_ff_equalization> change. In a running MATLAB session, also call
% <clear get_hrtf_equalization_filter> to reload the new file.
%
% OUTPUT (.mat file) containing:
%   ffL, ffR : [nTaps x 1] minimum-phase FIR filters (left, right ear) at fs
%   fs       : sampling frequency [Hz] (44100, FABIAN database)
%   f        : frequency vector of the magnitudes below [Hz]
%   H0L, H0R : magnitude of the frontal HRTF used as reference (linear)
%   GL, GR   : target magnitude of the filters (linear)
%   settings : design options (see <get_hrtf_ff_equalization>)
%   source   : provenance information
%
% Author: Gil Felix Greco, Braunschweig 11.11.2026
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to
%   develop this code under the author's direction. All codes were verified
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

% root folder of the framework (this script is in <root>/utilities)
root = fileparts( fileparts( mfilename('fullpath') ) );

%% design (default settings: method 'band', 50 Hz - 16 kHz, 20 dB boost limit)
FF = get_hrtf_ff_equalization();

ffL = FF.left;
ffR = FF.right;
fs  = FF.fs;
f   = FF.f;
H0L = FF.H0_left;   H0R = FF.H0_right;
GL  = FF.G_left;    GR  = FF.G_right;
settings = FF.settings;
source = sprintf( ['Free-field equalization filters of the FABIAN HRTFs (measured, SH ' ...
    'representation, HATO = 0): inverted frontal HRTF (az = 0, el = 0), per ear, method ''%s''. ' ...
    'Generated with utilities/generate_FABIAN_FF_equalization_mat.m and ' ...
    'auralization/hrtf_equalization/get_hrtf_ff_equalization.m on %s. ' ...
    'FABIAN database: Brinkmann et al., DOI 10.14279/depositonce-5718.5'], ...
    settings.method, datestr( now, 'yyyy-mm-dd' ) );

out_file = fullfile( root, 'third_party', 'FABIAN_HRTF_DATABASE_v4', ...
                     'FABIAN_FF_equalization_measured_HATO_0.mat' );
save( out_file, 'ffL', 'ffR', 'fs', 'f', 'H0L', 'H0R', 'GL', 'GR', 'settings', 'source', '-v7' );
fprintf( 'Saved %d-tap FF equalization filters (fs = %g Hz) to\n  %s\n', numel(ffL), fs, out_file );

%% check plot: frontal HRTF, FF filter and their product (should be 0 dB per band)
nfft = numel( f )*2 - 2;
[h0L, h0R] = AKhrirInterpolation( 0, 0, settings.HATO, 'measured_sh' );
EL = abs( fft( ffL, nfft ) );  EL = EL(1:end/2+1);
ER = abs( fft( ffR, nfft ) );  ER = ER(1:end/2+1);
HL = abs( fft( h0L, nfft ) );  HL = HL(1:end/2+1);
HR = abs( fft( h0R, nfft ) );  HR = HR(1:end/2+1);

figure( 'Name', 'FABIAN FF equalization' );
subplot( 2, 1, 1 );
semilogx( f, 20*log10(HL), f, 20*log10(EL), f, 20*log10(HL.*EL) );
grid on; xlim( [20 20000] ); ylim( [-30 30] );
ylabel( 'Magnitude (dB)' ); title( 'Left ear' );
legend( 'frontal HRTF |H_0|', 'FF filter', 'product (narrowband)', 'Location', 'southwest' );
subplot( 2, 1, 2 );
semilogx( f, 20*log10(HR), f, 20*log10(ER), f, 20*log10(HR.*ER) );
grid on; xlim( [20 20000] ); ylim( [-30 30] );
xlabel( 'Frequency (Hz)' ); ylabel( 'Magnitude (dB)' ); title( 'Right ear' );
