function hrtf = get_hrtf_database( database, subject )
% function hrtf = get_hrtf_database( database, subject )
%
% Returns the HRTF set used for the binaural rendering and for the HRTF equalizations
% (DF, FF). Without input arguments, the set is read from the .ini file:
%
%   hrtf_database = FABIAN (DEFAULT) or HUTUBS
%   hutubs_subject = 1 ... 96 (only used, and then required, with hrtf_database = HUTUBS)
%
% With input arguments, the set is defined directly (e.g. in post-processing scripts,
% without .ini file):
%   hrtf = get_hrtf_database( 'FABIAN' );
%   hrtf = get_hrtf_database( 'HUTUBS', 26 );
%
% Both databases provide HRTFs at the blocked ear-canal entrance, referenced to the sound
% pressure at the centre of the head with the head absent, at 44.1 kHz (256 samples), in
% the same head-centred coordinates (azimuth 0 = front, 90 = left ear; elevation 90 =
% above, 0 = front, -90 = below). Only the MEASURED HRTFs are used, with HATO = 0.
%
%   FABIAN: FABIAN head-and-torso simulator, measured in Oldenburg, SH order 35 (11950
%     directions). Shipped with the framework (third_party/FABIAN_HRTF_DATABASE_v4/).
%     HRIRs by <AKhrirInterpolation>. DF: inverted CTF provided with the database;
%     FF: precomputed filters (see <get_hrtf_equalization_filter>).
%     Brinkmann et al. (2017), J. Audio Eng. Soc. 65(10), 841-848,
%     DOI: 10.17743/jaes.2017.0033; data: DOI 10.14279/depositonce-5718.5
%
%   HUTUBS: 96 subjects (93 humans and the FABIAN head-and-torso simulator), measured at
%     TU Berlin, SH order 16 (440 directions). NOT shipped: download from
%     DOI 10.14279/depositonce-8487 and place the HRIR files in
%     third_party/HUTUBS_HRTF_DATABASE/ (only ppN_SHcoefficients_measured.mat is used).
%     HRIRs by <AKhrirDatabase> from the SH coefficients (loaded once per subject).
%     DF and FF filters are computed when first needed (<get_hrtf_df_equalization>,
%     <get_hrtf_ff_equalization>).
%     Subjects 1 and 96: FABIAN (repeated measurement), subjects 22 and 88: the same
%     human (repeated measurement). HUTUBS pp1/pp96 are NOT identical to
%     hrtf_database = FABIAN (other measurement system; differences below 1 dB up to
%     8 kHz and up to 3 dB at 20 kHz, see Brinkmann et al. 2019, Sect. 3.1).
%     Brinkmann et al. (2019), J. Audio Eng. Soc. 67(9), 705-718,
%     DOI: 10.17743/jaes.2019.0024; data: DOI 10.14279/depositonce-8487
%
% OUTPUT
%   hrtf : struct with fields
%     .database  : 'FABIAN' or 'HUTUBS'
%     .subject   : [] (FABIAN) or subject number (HUTUBS)
%     .name      : 'FABIAN' or 'HUTUBS pp<N>' (for printing)
%     .tag       : 'FABIAN' or 'HUTUBS_pp<N>' (key for the cached equalization filters)
%     .fs        : sampling frequency of the HRIRs (44100 Hz)
%     .file      : file with the HRTF data
%     .sh        : (HUTUBS only) SH coefficients, as loaded from .file
%     .reference : reference of the database (for printing)
%
% The HUTUBS data are cached (persistent). After replacing a database file in the same
% MATLAB session, run <clear get_hrtf_database>.
%
% Gil Felix Greco, Braunschweig 19.11.2026
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to develop this
%   code under the author's direction. All codes were verified by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global input_file
persistent cache

%% database and subject: from the input arguments or from the .ini file

if nargin == 0
    if isfield( input_file, 'hrtf_database' ) && ~isempty( strtrim( char( input_file.hrtf_database ) ) )
        database = input_file.hrtf_database;
    else
        database = 'FABIAN'; % default value
    end
    if isfield( input_file, 'hutubs_subject' )
        subject = input_file.hutubs_subject;
    else
        subject = [];
    end
elseif nargin < 2
    subject = [];
end

database = upper( strtrim( char( database ) ) );

switch database

    case 'FABIAN'

        hrtf.database  = 'FABIAN';
        hrtf.subject   = [];
        hrtf.name      = 'FABIAN';
        hrtf.tag       = 'FABIAN';
        hrtf.fs        = 44100;
        hrtf.file      = 'FABIAN_HRIR_measured_HATO_0.mat';
        hrtf.reference = ['FABIAN HRTF database (measured, HATO = 0), Brinkmann et al., ' ...
                          'J. Audio Eng. Soc. 65(10), 2017; DOI: 10.14279/depositonce-5718.5'];

    case 'HUTUBS'

        subject = il_check_subject( subject );

        hrtf.database  = 'HUTUBS';
        hrtf.subject   = subject;
        hrtf.name      = sprintf( 'HUTUBS pp%d', subject );
        hrtf.tag       = sprintf( 'HUTUBS_pp%d', subject );
        hrtf.fs        = 44100;
        hrtf.file      = sprintf( 'pp%d_SHcoefficients_measured.mat', subject );
        hrtf.reference = ['HUTUBS HRTF database (measured), Brinkmann et al., ' ...
                          'J. Audio Eng. Soc. 67(9), 2019; DOI: 10.14279/depositonce-8487'];

        % load the SH coefficients once per subject
        if isstruct( cache ) && isfield( cache, hrtf.tag )
            hrtf.sh = cache.(hrtf.tag);
        else
            if exist( hrtf.file, 'file' ) ~= 2
                error( 'auralization:hutubsNotFound', ...
                    ['<hrtf_database> = HUTUBS needs the file <%s>, which is not on the MATLAB path.\n' ...
                     'The HUTUBS database is not shipped with the framework. Download it from ' ...
                     'https://doi.org/10.14279/depositonce-8487 (file HRIRs.zip), extract it to ' ...
                     'third_party/HUTUBS_HRTF_DATABASE/ and run startup_auralization.m again. ' ...
                     'Only the files ppN_SHcoefficients_measured.mat are used.'], hrtf.file );
            end
            sh = load( hrtf.file );
            il_check_sh( sh, hrtf );
            cache.(hrtf.tag) = sh;
            hrtf.sh = sh;
        end

    otherwise
        error( 'auralization:invalidHrtfDatabase', ...
            '<hrtf_database> in the .ini file must be FABIAN or HUTUBS, found: %s', database );
end

end

%% inline function: check the HUTUBS subject number
function subject = il_check_subject( subject )

if isempty( subject )
    error( 'auralization:hutubsSubjectMissing', ...
        ['<hrtf_database> = HUTUBS needs <hutubs_subject> (1 to 96) in the .ini file. ' ...
         'Subjects 1 and 96 are the FABIAN head-and-torso simulator.'] );
end

if ischar( subject ) || isstring( subject )
    subject_str = strtrim( char( subject ) );
    subject     = str2double( regexprep( subject_str, '^pp', '', 'ignorecase' ) ); % accepts 26 or pp26
else
    subject_str = num2str( subject );
end

if ~isscalar( subject ) || ~isfinite( subject ) || subject ~= round( subject ) ...
        || subject < 1 || subject > 96
    error( 'auralization:invalidHutubsSubject', ...
        '<hutubs_subject> in the .ini file must be an integer from 1 to 96, found: %s', subject_str );
end

end

%% inline function: check the loaded SH coefficients
function il_check_sh( sh, hrtf )

ok = isfield( sh, 'HRIR' ) && isfield( sh.HRIR, 'CoeffLeft' ) && isfield( sh.HRIR, 'CoeffRight' ) ...
     && isfield( sh, 'TOA' ) && isfield( sh.HRIR, 'SamplingRate' );
if ~ok
    error( 'auralization:hutubsFileFormat', ...
        ['<%s> does not contain the expected fields (HRIR.CoeffLeft, HRIR.CoeffRight, ' ...
         'HRIR.SamplingRate, TOA). Use the MEASURED SH coefficients of the HUTUBS database.'], hrtf.file );
end

if double( sh.HRIR.SamplingRate ) ~= hrtf.fs
    error( 'auralization:hutubsSamplingFrequency', ...
        'Sampling frequency of <%s> is %g Hz, expected %g Hz.', ...
        hrtf.file, double( sh.HRIR.SamplingRate ), hrtf.fs );
end

end
