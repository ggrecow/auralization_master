function auralization_master(main_input_path_in, tag, input_file_path, results_path)
% function auralization_master(main_input_path_in, tag, input_file_path, results_path)
%
% Auralization framework using PANAM inputs to describe the sound source
% and flight profile. The codes performs the following:
%
% 1) inputs from PANAM are loaded, converted to matlab format and prepared
%     for auralization. The <main_input_path> needs to contain the following data files:
%
%   - Sound source: uses the 'auralization_input.dat' input file from PANAM 
%   - Flight profile: uses the 'geschw_hoehe_verlauf.dat' input file from PANAM
%
%   2) the sound source description is synthesised to a time-tomain sound
%       pressure signal, with a max. freq bounded to the half of the sampling rate 
%       (i.e. max. freq = fs/2) while the freq descritization, df, is bounded to
%       the inverse of the dt (i.e. df = 1/dt) which the PANAM predictions are provided Tonal and 
%       broadband noise components are synthesed separetly and them combined, providing 
%       sound files describing each sound source provided by PANAM. These are then combined 
%       to provided the synthesised sound source audio files of the engine, airframe and overall 
%       aircraft noise components.
%
%   3)  the atmospheric transfer function, which describes the freq. dependent attenuation 
%         imposed to the sound propagation by the atmosphere is calculated using the ART 
%         ray-tracing, which is already incorporated into this framework. For that purpose, an 
%         atmosphere is defined using the <input_file.m>, and the atmospheric transfer function 
%         is calculated for each source/receiver combination, as provided
%         by the  'geschw_hoehe_verlauf.dat'  (source positions) and the 
%         'auralization_input.dat' (receiver position). The sound
%         propagation effects considered are: geometrical spreading, air
%         attenuation, ground reflection (can be turned on/off). Doppler
%         effects needs to be already included in the inputs provided by
%         PANAM.         
%         
%   4) the calculated atmospheric transfer functions for each pair or
%       source/receiver positions are transformed into FIR filters and
%       applied to the synthesised sound source files to obtain the
%       auralized sound immission audio files (i.e. the sound files
%       describing the aircraft noise in the ground.)
%
%     Remarks:
%     - The number of receiver positions is automatically obtained from the
%     <auralization_input.dat> and the auralization is performed for all of
%     them. Therefore, don't use an input file with MANY observer positions
%     unless desired.
%
%      - Doppler effect needs to be already included in EMISSION input data
%      <'auralization_input.dat'>
%     
%      - It is preferred to use panam predictions 
%        with time-step 0.1 otherwise ground reflections are not smoothly synthesized. An interpolation  
%        of the flight trajectory used to compute the atmospheric transfer
%        functions could be used to avoid that, however it was choosen not to do that because
%        interpolating the flight trajectory is tricky as we have multiple dimensions 
%       (x,z,y,velocity, ...) depending on each other.
%
%   INPUTS 
%
%   main_input_path : string (mandatory)
%       main path pointing to the folder where the code can find the inputs
%       necessary to perform the auralizations from PANAM. The following inputs need to be 
%       found in this folder: 
%
%       -  <auralization_input.dat> input file from PANAM providing the sound
%           source description
%
%       - <geschw_hoehe_verlauf.dat> input file from PANAM providing the flight profile
%
%       The path can be given with or without a trailing file separator.
%
%    tag : string (optional)
%       provides the name of the case being analysed. All the figures and
%       sound files will be saved using this tag. If nothing is provided,
%       then a default string 'auralization_results' is used
%
%    input_file_path : string (optional)
%       contains all input parameters necessary to define the atmosphere and
%       perform to the auralizations. If the full path pointing to the <input_file.ini> is not provided
%       (or is empty), the code looks for an <*.ini> file inside the <main_input_path>.
%       WARNING: an error is thrown if no <*.ini> file or more than one <*.ini> file
%       is found inside <main_input_path>. In this case, provide <input_file_path> explicitly.
%
%   results_path : string (optional)
%   path where all the output data shall be saved. In this path, a folder called 
%   'auralization_results' will be created automatically, and the results for each receiver 
%   position founded in the input data fom PANAM will be saved within separate folders automatically. 
%   By default, results are saved within <main_input_path> if a <results_path> folder is not provided
%
% -------------------------------
% Author: Gil Felix Greco (ggrecow@gmail.com)
% Institution: Technische Universität Braunschweig 
%
% Date created: 11.04.2024
% Date last modified: 25.09.2026
% MATLAB version: 2024b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with input args

% remove trailing file separators (if any), so <main_input_path_in> can be given with or without them
main_input_path_in = regexprep( char(main_input_path_in), '[\\/]+$', '' );
main_input_path = [main_input_path_in filesep];

if nargin < 2 || isempty( tag )
    tag = 'auralization_results';
end

% if no <input_file_path> is provided, look for exactly one .ini file inside <main_input_path>
if nargin < 3 || isempty( input_file_path )
    input_file_path = find_ini_file( main_input_path_in );
elseif ~isfile( input_file_path )
    error( 'auralization:iniNotFound', ...
        'The <input_file_path> provided does not exist:\n  %s', input_file_path );
end

fprintf( 'Input file used: %s\n', input_file_path );

%% set global variables

% load input_file contaning key input parameteres
global input_file
input_file = ini2struct ( input_file_path ); 

% resolve relative paths given inside the .ini relative to the folder of the .ini
% (so they work regardless of MATLAB's current folder, also in the compiled .exe)
if isfield( input_file, 'sounding_filepath' ) && ~isempty( input_file.sounding_filepath )
    input_file.sounding_filepath = resolve_ini_path( input_file.sounding_filepath, fileparts(input_file_path) );
    fprintf( 'Atmospheric sounding file: %s\n', input_file.sounding_filepath );
end

% save figures in .fig?
global save_mat_fig
save_mat_fig = 0;

%%  show/save options    

% plot flight profile from panam
show_flight_profile = 1;  

% save_figs from flight profile: 0 (no); 1(yes)
flight_profile_save_fig = 1;  

% plot outputs from auralization
show_auralization = 1;     

% save figures from auralization
save_figs = 1;

%% Input data (sound source from PANAM, convert data to matlab)

% input file containing source description from PANAM
input_source = 'auralization_input.dat';
PATH_source = [main_input_path  input_source];  

% function : PANAM to SQAT data conversion + conventional noise metrics  
[source_data,...
 source_OASPL, source_OASPL_dBA,...
 source_SPECTROGRAM, source_SPECTROGRAM_dBA] = PANAM2struct( PATH_source ); 

%% create results folder

% number of receivers on the input data
nReceiver = size( source_data, 2 );

% By default, results are saved within <main_input_path> if a
% <results_path> folder is not provided
if nargin < 4 || isempty( results_path )
    results_path = main_input_path;
end

resultsFolder = create_auralization_results_folder( results_path, nReceiver );

%% Input data (flight profile from PANAM)

% name of the file containing the flight profile from PANAM
flight_profile_input = 'geschw_hoehe_verlauf.dat';

%  boolean flag : 0 (approach); 1 (departure); 2 (flyover) - only changes the x-axis label
flight_procedure = 2;   

% function : plot/show/save flight profile
flight_profile = get_flight_profile( [main_input_path flight_profile_input],... % data path
                                                                                 show_flight_profile,... % plot flight profile ? 0 (no); 1(yes)
                                                                                     flight_procedure,... % proc: 0 (approach); 1 (departure); 2 (flyover) - only changes the x-axis label
                                                                            flight_profile_save_fig,... % save_figs: 0 (no); 1(yes)
                                                                                                          tag,... % tag: string with name tag to save figure
                                                                    resultsFolder.mainFolder);   % figuresdir: string with dir path to save figs
clc;

clear  show_flight_profile flight_profile_input flight_procedure flight_profile_save_fig;

%% main auralization code for each receiver position

% time after and before the max SPL on the input data
% trimTime = 20;  % seconds
trimTime = str2double( input_file.trim_time );

% number of smoothing iteration during the conversion from toc to narrowband (can be ZERO)
% smoothings = 0;  

% By default, smoothings_emission_based = 0, so it doesnt even need to be
% provided in the <input_file>
if isfield( input_file, 'smoothings_emission_based' )
    smoothings_emission_based = str2double( input_file.smoothings_emission_based );  % get <smoothings_emission_based> from <input_file>
else
    smoothings_emission_based = 0; % not really necessary when auralizing using ray tracing propagation
end

for i = 1:nReceiver
    %----------------------------------------------------------------------------------------
    % deals with EMISSION data   
    %----------------------------------------------------------------------------------------

    if save_figs == 1
        tag_auralization = [resultsFolder.receiverFolder{i}  tag '_Receiver_' sprintf('%01d',i)];
    else
        tag_auralization = [];
    end
    
    % function : trim data before SQ metics calculation / auralization 
    [source_data_trimmed,...
     source_SPECTROGRAM_trimmed,...
     source_SPECTROGRAM_dBA_trimmed,...
     flight_profile_trimmed] = prepare_input_SQ( source_data(:,i),...                                % input_1
                                                                            source_SPECTROGRAM(1,i),...          % input_2
                                                                            source_SPECTROGRAM_dBA(1,i),... % input_3
                                                                            flight_profile,...                                      % flight profile
                                                                            trimTime, tag_auralization);
   
     % call main auralization function (i.e. sinthesize source data in time-domain)        
     input_type = 'emission';

     % auralization based on engine/airframe noise sources
     MASTER_auralization_EngineAirframe( source_data_trimmed, ...
                                                                                                             flight_profile_trimmed, ...
                                                                                                             smoothings_emission_based, ...
                                                                                                             input_type, ...
                                                                                                             tag_auralization, ...
                                                                                                             show_auralization );

     close all;

 end

    function ini_path = find_ini_file(folder)
        % Returns the full path of the only .ini file inside <folder>.
        % Throws an error if no .ini file or more than one .ini file is found.
        ini_candidates = dir( fullfile(folder, '*.ini') );
        switch numel( ini_candidates )
            case 0
                error( 'auralization:noIni', ...
                    ['No .ini file found inside <main_input_path>:\n  %s\n' ...
                     'Please provide it as third input: auralization_master(main_input_path, tag, input_file_path)'], ...
                    folder );
            case 1
                ini_path = fullfile( ini_candidates.folder, ini_candidates.name );
                fprintf( 'No <input_file_path> provided. Using the only .ini file found inside <main_input_path>.\n' );
            otherwise
                error( 'auralization:ambiguousIni', ...
                    ['More than one .ini file found inside <main_input_path>:\n  %s\n\n  %s\n\n' ...
                     'Please specify which one to use as third input: auralization_master(main_input_path, tag, input_file_path)'], ...
                    folder, strjoin( {ini_candidates.name}, [newline '  '] ) );
        end
    end % end function <find_ini_file>

    function p = resolve_ini_path(p, ini_folder)
        % Converts a path given inside the .ini into a full path.
        % Absolute paths are kept; relative paths are resolved relative to the .ini folder.
        p = strtrim( erase(p, {'"', ''''}) );                             % remove quotes, if any
        is_absolute = ~isempty( regexp(p, '^([a-zA-Z]:[\\/]|[\\/])', 'once') );
        if ~is_absolute
            p = fullfile( ini_folder, p );
        end
    end % end function <resolve_ini_path>

    function Result = ini2struct(FileName)
        %==========================================================================
        %  Author: Andriy Nych ( nych.andriy@gmail.com )
        % Version:        733341.4155741782200
        %==========================================================================
        %
        % INI = ini2struct(FileName)
        %
        % This function parses INI file FileName and returns it as a structure with
        % section names and keys as fields.
        %
        % Sections from INI file are returned as fields of INI structure.
        % Each fiels (section of INI file) in turn is structure.
        % It's fields are variables from the corresponding section of the INI file.
        %
        % If INI file contains "oprhan" variables at the beginning, they will be
        % added as fields to INI structure.
        %
        % Lines starting with ';' and '#' are ignored (comments).
        %
        % See example below for more information.
        %
        % Usually, INI files allow to put spaces and numbers in section names
        % without restrictions as long as section name is between '[' and ']'.
        % It makes people crazy to convert them to valid Matlab variables.
        % For this purpose Matlab provides GENVARNAME function, which does
        %  "Construct a valid MATLAB variable name from a given candidate".
        % See 'help genvarname' for more information.
        %
        % The INI2STRUCT function uses the GENVARNAME to convert strange INI
        % file string into valid Matlab field names.
        %
        % [ test.ini ]~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %
        %     SectionlessVar1=Oops
        %     SectionlessVar2=I did it again ;o)
        %     [Application]
        %     Title = Cool program
        %     LastDir = c:\Far\Far\Away
        %     NumberOFSections = 2
        %     [1st section]
        %     param1 = val1
        %     Param 2 = Val 2
        %     [Section #2]
        %     param1 = val1
        %     Param 2 = Val 2
        %
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %
        % The function converts this INI file it to the following structure:
        %
        % [ MatLab session (R2006b) ]~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %  >> INI = ini2struct('test.ini');
        %  >> disp(INI)
        %         sectionlessvar1: 'Oops'
        %         sectionlessvar2: 'I did it again ;o)'
        %             application: [1x1 struct]
        %             x1stSection: [1x1 struct]
        %            section0x232: [1x1 struct]
        %
        %  >> disp(INI.application)
        %                    title: 'Cool program'
        %                  lastdir: 'c:\Far\Far\Away'
        %         numberofsections: '2'
        %
        %  >> disp(INI.x1stSection)
        %         param1: 'val1'
        %         param2: 'Val 2'
        %
        %  >> disp(INI.section0x232)
        %         param1: 'val1'
        %         param2: 'Val 2'
        %
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %
        % NOTE.
        % WhatToDoWithMyVeryCoolSectionAndVariableNamesInIniFileMyVeryCoolProgramWrites?
        % GENVARNAME also does the following:
        %   "Any string that exceeds NAMELENGTHMAX is truncated". (doc genvarname)
        % Period.
        %
        % =========================================================================
        Result = [];                            % we have to return something
        CurrMainField = '';                     % it will be used later
        f = fopen(FileName,'r');                % open file
        if f < 0
            error( 'auralization:iniNotReadable', ...
                'Could not open the .ini file:\n  %s', FileName );
        end
        while ~feof(f)                          % and read until it ends
            s = strtrim(fgetl(f));              % Remove any leading/trailing spaces
            if isempty(s)
                continue;
            end;
            if (s(1)==';')                      % ';' start comment lines
                continue;
            end;
            if (s(1)=='#')                      % '#' start comment lines
                continue;
            end;
            if ( s(1)=='[' ) && (s(end)==']' )
                % We found section
                CurrMainField = genvarname(lower(s(2:end-1)));
                Result.(CurrMainField) = [];    % Create field in Result
            else
                % ??? This is not a section start
                [par,val] = strtok(s, '=');
                val = CleanValue(val);
                if ~isempty(CurrMainField)
                    % But we found section before and have to fill it
                    Result.(CurrMainField).(lower(genvarname(par))) = val;
                else
                    % No sections found before. Orphan value
                    Result.(lower(genvarname(par))) = val;
                end
            end
        end
        fclose(f);
        return;
        function res = CleanValue(s)
            res = strtrim(s);
            if strcmpi(res(1),'=')
                res(1)=[];
            end
            res = strtrim(res);
            return;
        end
    end % end function <ini2struct>

end