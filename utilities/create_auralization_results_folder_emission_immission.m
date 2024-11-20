function resultsFolder = create_auralization_results_folder_emission_immission(main_input_path, nReceiver)
%function resultsFolder = create_auralization_results_folder_emission_immission(main_input_path, nReceiver)
%
% This function creates the folders to store the results from the
% auralization procedure (both emission-based and Immission-based). It creates the folder inside the
% <main_input_path> (the input folder path is conventioned/assumed to end
% with \).
%
% The following folders are created:
%
%   - auralization_results : main results folder to include all results
%     contained within the input file 
%
%        -- Emission_based : main results folder to include all results of
%           auralization procedure conducted using the emission data from panam as input
%
%           --- Receiver : subsub folder created inside <figs> folder to
%                store the figures and auralized audio files in subfolders for each receiver position
%
%   - auralization_results : main results folder to include all results
%     contained within the input file 
%
%        -- Immission_based : main results folder to include all results of
%           auralization procedure conducted using the immision data from panam as input
%
%           --- Receiver : subsub folder created inside <figs> folder to
%                store the figures and auralized audio files in subfolders for each receiver position
%
%   OUTPUT : struct
%   resultsFolder contains all folder paths within it
%
%   resultsFolder.mainFolder
%   resultsFolder.emissionBased
%   resultsFolder.receiverFolder_emissionBased (one for each receiver position)
%   resultsFolder.immissionBased
%   resultsFolder.receiverFolder_imissionBased (one for each receiver position)
%
% Author: Gil Felix Greco, Braunschweig 24.04.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% create main <results_auralization> folder within <main_input_path>

resultsFolder.mainFolder = [ main_input_path 'results_auralization\' ];

if ~exist( resultsFolder.mainFolder,'dir' )
    mkdir( resultsFolder.mainFolder );
end

%% create subsubfolder <Emission_based> inside the folder <results_auralization>

    resultsFolder.emissionBased = [ resultsFolder.mainFolder 'Emission_based\' ];

    if ~exist( resultsFolder.emissionBased , 'dir' )
        mkdir( resultsFolder.emissionBased );
    end

%% create subsubfolder <Receiver_nReceiver> inside the folder <Emission_based>>

for i = 1:nReceiver

    resultsFolder.receiverFolder_emissionBased{i} = [ resultsFolder.emissionBased 'Receiver_' sprintf('%01d',i) '\' ];

    if ~exist( resultsFolder.receiverFolder_emissionBased{i} , 'dir' )
        mkdir( resultsFolder.receiverFolder_emissionBased{i} );
    end

end

%% create subsubfolder <Immission_based> inside the folder <results_auralization>

    resultsFolder.immissionBased = [ resultsFolder.mainFolder 'Immission_based'  '\' ];

    if ~exist( resultsFolder.immissionBased , 'dir' )
        mkdir( resultsFolder.immissionBased );
    end

%% create subsubfolder <Receiver_nReceiver> inside the folder <Immission_based>>

for i = 1:nReceiver

    resultsFolder.receiverFolder_immissionBased{i} = [ resultsFolder.immissionBased 'Receiver_' sprintf('%01d',i) '\' ];

    if ~exist( resultsFolder.receiverFolder_immissionBased{i} , 'dir' )
        mkdir( resultsFolder.receiverFolder_immissionBased{i} );
    end

end

end