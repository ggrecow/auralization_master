function resultsFolder = create_auralization_results_folder(main_input_path, nReceiver)
%function resultsFolder = create_auralization_results_folder(main_input_path, nReceiver)
%
% This function creates the folders to store the results from the
% auralization procedure. It creates the folders inside the
% <main_input_path>, which can be given with or without a trailing file
% separator. Paths are built with <fullfile>, so the function works on
% Windows, Linux and macOS.
%
% The following folders are created:
%
%   - results_auralization : main results folder to include all results
%     contained within the input file
%
%      -- Receiver_<n> : sub folder created inside <results_auralization> to
%          store the figures and auralized audio files of each receiver position
%
%   OUTPUT : struct
%   resultsFolder contains all folder paths within it
%
%   resultsFolder.mainFolder
%   resultsFolder.receiverFolder (one for each receiver position)
%
%   All returned paths END WITH A FILE SEPARATOR, because the calling functions
%   append file names directly to them (e.g. [resultsFolder.receiverFolder{i} tag ...]).
%
% Author: Gil Felix Greco, Braunschweig 11.04.2024
% Date last modified: 27.10.2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% create main <results_auralization> folder within <main_input_path>

resultsFolder.mainFolder = [ fullfile( char(main_input_path), 'results_auralization' ) filesep ];

if ~exist( resultsFolder.mainFolder,'dir' )
    mkdir( resultsFolder.mainFolder );
end

%% create subsubfolder <Receiver_nReceiver> inside the folder <results_auralization>

resultsFolder.receiverFolder = cell( 1, nReceiver );

for i = 1:nReceiver

    resultsFolder.receiverFolder{i} = [ fullfile( resultsFolder.mainFolder, sprintf('Receiver_%01d', i) ) filesep ];

    if ~exist( resultsFolder.receiverFolder{i} , 'dir' )
        mkdir( resultsFolder.receiverFolder{i} );
    end

end

end