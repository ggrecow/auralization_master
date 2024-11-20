function resultsFolder = create_auralization_results_folder(main_input_path, nReceiver)
%function resultsFolder = create_auralization_results_folder(main_input_path, nReceiver)
%
% This function creates the folders to store the results from the
% auralization procedure. It creates the folder inside the
% <main_input_path> (the input folder path is conventioned/assumed to end
% with \).
%
% The following folders are created:
%
%   - auralization_results : main results folder to include all results
%     contained within the input file 
%
%      -- Receiver : sub folder created inside <figs> folder to
%          store the figures and auralized audio files in subfolders for each receiver position
%
%   OUTPUT : struct
%   resultsFolder contains all folder paths within it
%
%   resultsFolder.mainFolder
%   resultsFolder.receiverFolder (one for each receiver position)
%
% Author: Gil Felix Greco, Braunschweig 11.04.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% create main <results_auralization> folder within <main_input_path>

resultsFolder.mainFolder = [ main_input_path 'results_auralization\' ];

if ~exist( resultsFolder.mainFolder,'dir' )
    mkdir( resultsFolder.mainFolder );
end

%% create subsubfolder <Receiver_nReceiver> inside the folder <results_auralization> 
 
for i = 1:nReceiver
    
    resultsFolder.receiverFolder{i} = [ resultsFolder.mainFolder 'Receiver_' sprintf('%01d',i) '\' ];

    if ~exist( resultsFolder.receiverFolder{i} , 'dir' )
        mkdir( resultsFolder.receiverFolder{i} );
    end

end
 
end