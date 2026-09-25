% convert_FABIAN_CTF_sofa2mat.m
%
% MAINTAINER SCRIPT (not needed by users): converts the inverted, 3rd-octave
% smoothed, minimum-phase common transfer function (CTF) of the FABIAN HRTF
% database from SOFA to .mat. The resulting .mat file is shipped with the
% framework, so that get_FIR.m can apply diffuse-field equalization without
% requiring the full database or the SOFA API (e.g. in the compiled .exe).
% This script is kept for documentation/reproducibility of that file.
%
% The CTF is provided by the database authors (Brinkmann et al.) and is the
% same filter used by AKhrirInterpolation(..., 'dir').
%
% REQUIREMENTS (only for running this script)
%   - SOFA API for Matlab/Octave (SOFAstart must have been called)
%   - FABIAN_CTF_measured_inverted_smoothed.sofa from the FABIAN database,
%     doi.org/10.14279/depositonce-5718
%
% OUTPUT
%   third_party/FABIAN_HRTF_DATABASE_v4/FABIAN_CTF_measured_inverted_smoothed.mat
%   containing:
%     ctf    : [N x 1] inverted CTF impulse response
%     fs     : sampling frequency [Hz]
%     source : provenance information
%
% Gil Felix Greco, Braunschweig 25.09.2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

% root folder of the local copy of the FABIAN database (adapt if needed)
db_root = 'C:\Users\greco\Desktop\FABIAN_HRTF_DATABASE_V1';

% find the CTF file anywhere inside the database folder
sofa_name = 'FABIAN_CTF_measured_inverted_smoothed.sofa';
hits = dir( fullfile( db_root, '**', sofa_name ) );

if isempty( hits )
    error( 'No <%s> found inside <%s>.', sofa_name, db_root );
elseif numel( hits ) > 1
    error( 'More than one <%s> found inside <%s>:\n%s', sofa_name, db_root, ...
           strjoin( fullfile( {hits.folder}, {hits.name} ), '\n' ) );
end

sofa_file = fullfile( hits(1).folder, hits(1).name );
fprintf( 'Using CTF file:\n  %s\n', sofa_file );

% output folder (next to the HRTF data used by AKhrirInterpolation)
this_folder = fileparts( mfilename('fullpath') );
out_file    = fullfile( this_folder, '..', 'third_party', 'FABIAN_HRTF_DATABASE_v4', ...
                        'FABIAN_CTF_measured_inverted_smoothed.mat' );

%% load SOFA object
Obj = SOFAload( sofa_file );

IR = squeeze( Obj.Data.IR );   % SOFA convention: [M x R x N]
if ~isvector( IR )
    error( 'Expected a single CTF impulse response, got an array of size %s.', mat2str(size(IR)) );
end

ctf = IR(:);
fs  = Obj.Data.SamplingRate;

source = sprintf( ['FABIAN HRTF database, %s, converted from SOFA on %s. ' ...
                   'Title: %s. License: %s'], ...
                   sofa_name, datestr(now, 'yyyy-mm-dd'), ...
                   Obj.GLOBAL_Title, Obj.GLOBAL_License );

%% save
save( out_file, 'ctf', 'fs', 'source' );
fprintf( 'Saved %d-tap CTF (fs = %g Hz) to\n  %s\n', numel(ctf), fs, out_file );

%% check plot: magnitude response of the inverted CTF
nfft = 2^nextpow2( max( 4096, numel(ctf) ) );
H    = fft( ctf, nfft );
f    = (0:nfft/2)' * fs / nfft;

figure;
semilogx( f(2:end), 20*log10( abs( H(2:nfft/2+1) ) ) );
grid on; xlim([20 fs/2]);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('FABIAN inverted CTF (diffuse-field equalization filter)');