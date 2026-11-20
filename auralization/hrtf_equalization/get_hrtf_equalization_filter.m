function EQ = get_hrtf_equalization_filter( eqType, hrtf )
% function EQ = get_hrtf_equalization_filter( eqType, hrtf )
%
% Returns the HRTF equalization filters of an HRTF set (FABIAN or HUTUBS, see
% <get_hrtf_database>), at the database sampling frequency (44.1 kHz), i.e. the
% filters that would be applied to the HRIRs before their resampling to the
% framework sampling frequency. Used by <convert_binaural_equalization>.
%
%   eqType = 'none' -> unit impulse (measured HRTFs, blocked ear-canal entrance)
%   eqType = 'DF'   -> diffuse-field equalization: inverted, 3rd-octave smoothed,
%                      minimum-phase common transfer function (CTF). Same filter
%                      for both ears.
%                      FABIAN: CTF provided with the FABIAN database,
%                              file FABIAN_CTF_measured_inverted_smoothed.mat
%                      HUTUBS: no CTF is provided with the database; designed from
%                              the SH coefficients with <get_hrtf_df_equalization>
%                              (same method; reproduces the FABIAN CTF within 0.2 dB
%                              per third-octave band, 50 Hz - 16 kHz)
%   eqType = 'FF'   -> free-field equalization: inverted frontal HRTF
%                      (az = 0, el = 0, HATO = 0), one filter per ear, see
%                      <get_hrtf_ff_equalization>.
%                      FABIAN: file FABIAN_FF_equalization_measured_HATO_0.mat
%                              (designed on the fly, with a warning, if missing)
%                      HUTUBS: designed on the fly
%
% The FABIAN .mat files are shipped in third_party/FABIAN_HRTF_DATABASE_v4/, which
% must be on the MATLAB path.
%
% INPUTS
%   eqType : 'none', 'DF' or 'FF'
%   hrtf   : struct from <get_hrtf_database> (DEFAULT: get_hrtf_database(), i.e. the
%            set of the .ini file, or FABIAN if no .ini file was read). In
%            post-processing scripts, pass the set explicitly, e.g.
%            get_hrtf_database( 'HUTUBS', 26 ).
%
% OUTPUT
%   EQ : struct with fields
%     .type          : 'none', 'DF' or 'FF'
%     .hrtf          : name of the HRTF set (e.g. 'FABIAN', 'HUTUBS pp26')
%     .left, .right  : [nTaps x 1] FIR filters at EQ.fs
%     .fs            : sampling frequency of the filters (44100 Hz)
%     .source        : file the filters were loaded from (or a description)
%
% Results are cached (persistent) per HRTF set and equalization type, since the
% filters do not depend on the auralization case (computing the HUTUBS filters
% takes well below 1 s). After regenerating a .mat file in the same MATLAB
% session, run <clear get_hrtf_equalization_filter> to reload it.
%
% Author: Gil Felix Greco, Braunschweig 11.11.2026
% Gil Felix Greco, Braunschweig 19.11.2026 - HRTF set as input (FABIAN or HUTUBS);
%   cache per HRTF set and equalization type
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to
%   develop this code under the author's direction. All codes were verified
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent cache

if nargin < 2 || isempty( hrtf )
    hrtf = get_hrtf_database();
end

key = [ hrtf.tag '_' eqType ];   % e.g. FABIAN_DF, HUTUBS_pp26_FF
if isstruct( cache ) && isfield( cache, key )
    EQ = cache.(key);
    return
end

EQ.type = eqType;
EQ.hrtf = hrtf.name;
EQ.fs   = hrtf.fs;

switch eqType

    case 'none'
        EQ.left   = 1;
        EQ.right  = 1;
        EQ.source = 'no equalization (unit impulse)';

    case 'DF'
        switch hrtf.database
            case 'FABIAN'
                ctf_file = 'FABIAN_CTF_measured_inverted_smoothed.mat';
                if ~exist( ctf_file, 'file' )
                    error( 'auralization:ctfNotFound', ...
                        ['Diffuse-field equalization requires <%s>, which is shipped in ' ...
                        'third_party/FABIAN_HRTF_DATABASE_v4/. Make sure this folder is on the MATLAB path.'], ctf_file );
                end
                CTF = load( ctf_file, 'ctf', 'fs' );
                il_check_fs( CTF.fs, hrtf.fs, ctf_file );
                EQ.left   = CTF.ctf(:);
                EQ.right  = CTF.ctf(:);
                EQ.source = ctf_file;
            otherwise
                DF = get_hrtf_df_equalization( hrtf );
                EQ.left   = DF.left(:);
                EQ.right  = DF.right(:);
                EQ.source = sprintf( ['inverted CTF of %s, designed from the SH coefficients ' ...
                                      '(get_hrtf_df_equalization)'], hrtf.name );
        end

    case 'FF'
        switch hrtf.database
            case 'FABIAN'
                ff_file = 'FABIAN_FF_equalization_measured_HATO_0.mat';
                if exist( ff_file, 'file' )
                    FF = load( ff_file, 'ffL', 'ffR', 'fs' );
                    il_check_fs( FF.fs, hrtf.fs, ff_file );
                    EQ.left   = FF.ffL(:);
                    EQ.right  = FF.ffR(:);
                    EQ.source = ff_file;
                else
                    warning( 'auralization:ffFileNotFound', ...
                        ['<%s> not found on the MATLAB path: the free-field equalization filters ' ...
                         'are designed on the fly (default settings of <get_hrtf_ff_equalization>).'], ff_file );
                    FF = get_hrtf_ff_equalization( struct( 'hrtf', hrtf ) );
                    EQ.left   = FF.left(:);
                    EQ.right  = FF.right(:);
                    EQ.source = 'get_hrtf_ff_equalization (designed on the fly)';
                end
            otherwise
                FF = get_hrtf_ff_equalization( struct( 'hrtf', hrtf ) );
                EQ.left   = FF.left(:);
                EQ.right  = FF.right(:);
                EQ.source = sprintf( ['inverted frontal HRTF of %s, designed on the fly ' ...
                                      '(get_hrtf_ff_equalization)'], hrtf.name );
        end

    otherwise
        error( 'auralization:invalidHrtfEqualization', ...
            'Unknown HRTF equalization type ''%s'' (use ''none'', ''DF'' or ''FF'').', eqType );
end

cache.(key) = EQ;

end

%% inline function: check sampling frequency of a filter file
function il_check_fs( fs_file, fs_hrir_db, fileName )
if fs_file ~= fs_hrir_db
    error( 'auralization:eqSamplingFrequency', ...
        'Sampling frequency of <%s> (%g Hz) differs from the HRTF database (%g Hz).', ...
        fileName, fs_file, fs_hrir_db );
end
end
