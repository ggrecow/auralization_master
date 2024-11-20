function sim_result_files = run_dynamic_scenario(obj,...
            source_positions, receiver_positions,...
            output_folder, force)
%METHOD1 Summary of this method goes here
%   Detailed explanation goes here

arguments
    obj (1,1) itaPigeonProject
    source_positions (:,3) double
    receiver_positions (:,3) double
    output_folder {mustBeTextScalar}
    force (1,1) logical = false
end

% Input checks
nSourcePos = size(source_positions, 1);
nReceiverPos = size(receiver_positions, 1);
nFrames = max([nSourcePos, nReceiverPos]);
assert( nSourcePos == 1 || nReceiverPos == 1 ||  nReceiverPos == nSourcePos,...
    'Number of source / receiver positions must be either one or match the respective other number.');

%Duplicate positions if required
if nReceiverPos == 1
    receiver_positions = repmat(receiver_positions, nSourcePos, 1);
end
if nSourcePos == 1
    source_positions = repmat(source_positions, nReceiverPos, 1);
end

% Output files
sim_result_files = GetSimFileList(output_folder, nFrames);
if ~exist(output_folder, 'dir')
     mkdir(output_folder); 
end

% Run urban propagation
if force || ~isfile( sim_result_files{end} ) || ~ismember('urbanPaths', who('-file', sim_result_files{end}) )
    urbanResultFileBackup = obj.result_file_path;
    urbanConfFileBackup = obj.config_file_path;
    urbanStatsFileBackup = obj.runtime_stats_path;
    run_quiet = obj.run_quiet;

    obj.result_file_path = fullfile(output_folder, 'pigeon_urban_paths.json');
    obj.config_file_path = fullfile(output_folder, 'pigeon_project.ini');
    obj.runtime_stats_path = fullfile(output_folder, 'pigeon_project_stats.json');
    obj.run_quiet = true;

    wb = itaWaitbar(nFrames,'Running urban propagation...');
    for idFrame = 1:nFrames
        urbanPaths = obj.run( source_positions(idFrame, :), receiver_positions(idFrame, :) );
        save(sim_result_files{idFrame}, 'urbanPaths'); %, '-append'
        if isvalid(wb); wb.inc(); end
    end
    if isvalid(wb); wb.close(); end
    
    obj.result_file_path = urbanResultFileBackup;
    obj.config_file_path = urbanConfFileBackup;
    obj.runtime_stats_path = urbanStatsFileBackup;
    obj.run_quiet = run_quiet;
end
end

function [frame_files] = GetSimFileList(outputFolder, nFrames)

    frame_files = cell(1, nFrames);
    for idFrame = 1:nFrames  
        frame_files{idFrame} = fullfile( outputFolder, ['pigeon_simulation_frame_',num2str(idFrame),'.mat'] );
    end
end