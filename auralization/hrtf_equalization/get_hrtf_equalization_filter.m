function EQ = get_hrtf_equalization_filter( eqType )
% function EQ = get_hrtf_equalization_filter( eqType )
%
% Returns the HRTF equalization filters of the FABIAN database, at the
% database sampling frequency (44.1 kHz), i.e. the filters that would be
% applied to the HRIRs before their resampling to the framework sampling
% frequency. Used by <convert_binaural_equalization>.
%
%   eqType = 'none' -> unit impulse (measured HRTFs, blocked ear-canal entrance)
%   eqType = 'DF'   -> diffuse-field equalization: inverted, 3rd-octave smoothed,
%                      minimum-phase common transfer function (CTF) provided with
%                      the FABIAN database. Same filter for both ears.
%                      File: FABIAN_CTF_measured_inverted_smoothed.mat
%   eqType = 'FF'   -> free-field equalization: inverted frontal HRTF
%                      (az = 0, el = 0, HATO = 0), one filter per ear, see
%                      <get_hrtf_ff_equalization>.
%                      File: FABIAN_FF_equalization_measured_HATO_0.mat
%
% Both .mat files are shipped in third_party/FABIAN_HRTF_DATABASE_v4/, which
% must be on the MATLAB path. If the FF file is missing, the filters are
% designed on the fly with the default settings of <get_hrtf_ff_equalization>
% (with a warning), which needs AKtools and the FABIAN SH data.
%
% OUTPUT
%   EQ : struct with fields
%     .type          : 'none', 'DF' or 'FF'
%     .left, .right  : [nTaps x 1] FIR filters at EQ.fs
%     .fs            : sampling frequency of the filters (44100 Hz)
%     .source        : file the filters were loaded from (or a description)
%
% Results are cached (persistent), since the filters do not depend on the
% auralization case. After regenerating a .mat file in the same MATLAB
% session, run <clear get_hrtf_equalization_filter> to reload it.
%
% Author: Gil Felix Greco, Braunschweig 11.11.2026
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to
%   develop this code under the author's direction. All codes were verified
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent cache
if isstruct( cache ) && isfield( cache, eqType )
    EQ = cache.(eqType);
    return
end

fs_hrir_db = 44100; % sampling frequency of the FABIAN HRTF database

EQ.type = eqType;
EQ.fs   = fs_hrir_db;

switch eqType

    case 'none'
        EQ.left   = 1;
        EQ.right  = 1;
        EQ.source = 'no equalization (unit impulse)';

    case 'DF'
        ctf_file = 'FABIAN_CTF_measured_inverted_smoothed.mat';
        if ~exist( ctf_file, 'file' )
            error( 'auralization:ctfNotFound', ...
                ['Diffuse-field equalization requires <%s>, which is shipped in ' ...
                'third_party/FABIAN_HRTF_DATABASE_v4/. Make sure this folder is on the MATLAB path.'], ctf_file );
        end
        CTF = load( ctf_file, 'ctf', 'fs' );
        il_check_fs( CTF.fs, fs_hrir_db, ctf_file );
        EQ.left   = CTF.ctf(:);
        EQ.right  = CTF.ctf(:);
        EQ.source = ctf_file;

    case 'FF'
        ff_file = 'FABIAN_FF_equalization_measured_HATO_0.mat';
        if exist( ff_file, 'file' )
            FF = load( ff_file, 'ffL', 'ffR', 'fs' );
            il_check_fs( FF.fs, fs_hrir_db, ff_file );
            EQ.left   = FF.ffL(:);
            EQ.right  = FF.ffR(:);
            EQ.source = ff_file;
        else
            warning( 'auralization:ffFileNotFound', ...
                ['<%s> not found on the MATLAB path: the free-field equalization filters ' ...
                 'are designed on the fly (default settings of <get_hrtf_ff_equalization>).'], ff_file );
            FF = get_hrtf_ff_equalization( struct( 'fs_hrir', fs_hrir_db ) );
            EQ.left   = FF.left(:);
            EQ.right  = FF.right(:);
            EQ.source = 'get_hrtf_ff_equalization (designed on the fly)';
        end

    otherwise
        error( 'auralization:invalidHrtfEqualization', ...
            'Unknown HRTF equalization type ''%s'' (use ''none'', ''DF'' or ''FF'').', eqType );
end

cache.(eqType) = EQ;

end

%% inline function: sampling frequency check
function il_check_fs( fs_file, fs_hrir_db, fileName )
if fs_file ~= fs_hrir_db
    error( 'auralization:eqSamplingFrequency', ...
        'Sampling frequency of <%s> (%g Hz) differs from the HRTF database (%g Hz).', ...
        fileName, fs_file, fs_hrir_db );
end
end
