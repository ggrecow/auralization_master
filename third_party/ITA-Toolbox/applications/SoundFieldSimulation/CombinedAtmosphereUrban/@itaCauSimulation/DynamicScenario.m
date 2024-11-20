function sim_result_files = DynamicScenario(obj, sourcePositions, receiverPositions, updateRate, outputFolder, force)
%DYNAMICSCENARIO Simulates a dynamic scenario using a set of source and
%receiver positions. The result for each frame are stored in a .mat files
%and a cell with respective filenames is returned.
%   Inputs (default):
%   sourcePositions:    Nx3 matrix with source positions [m]
%   receiverPositions:  Nx3 matrix with receiver positions [m]
%   updateRate:         Update rate of positions (= blocksize / samplingrate)
%   outputFolder:       Folder used to store the simulation results
%   force (false):      If true, simulation is carried out even if files already exist [boolean]
%
%   It is also possible to use a static source / receiver by only handing a
%   1x3 vector to this function.

if nargin < 6; force = false; end

%% Input checks
assert( size(sourcePositions,2) == 3, 'Source Positions must be defined as Nx3-matrix')
assert( size(receiverPositions,2) == 3, 'Receiver Positions must be defined as Nx3-matrix')

nSourcePos = size(sourcePositions, 1);
nReceiverPos = size(receiverPositions, 1);
nFrames = max([nSourcePos, nReceiverPos]);
assert( nSourcePos == 1 || nReceiverPos == 1 ||  nReceiverPos == nSourcePos,...
    'Number of source / receiver positions must be either one or match the respective other number.');

%Duplicate positions if required
if nReceiverPos == 1
    receiverPositions = repmat(receiverPositions, nSourcePos, 1);
end
if nSourcePos == 1
    sourcePositions = repmat(sourcePositions, nReceiverPos, 1);
end

%% Output files
[sim_result_files, delay_file] = GetFileList(outputFolder, nFrames);
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end

%% Atmospheric Ray Tracing
if force || ~isfile( sim_result_files{end} ) || ~ismember('atmoPath', who('-file', sim_result_files{end}) )
    
    wb = itaWaitbar(nFrames,'Running Ray Tracing...');
    for idFrame = 1:nFrames  
        [atmoPath, receiverMissed] = obj.RunRayTracing(sourcePositions(idFrame,:),receiverPositions(idFrame, :));
        if receiverMissed
            warning('Eigenrays could not be determined. Using previous path.')
            load(sim_result_files{idFrame-1}, 'atmoPath');
        end
        save(sim_result_files{idFrame},'atmoPath');
        if isvalid(wb); wb.inc(); end
    end
    
    if isvalid(wb); wb.close(); end
end

%% Extract propagation delay and apply low pass
%Check if stored data matches number of frames
if isfile( delay_file )
    load(delay_file, 'propDelay');
    nFramesDelay = numel(propDelay);
end

if force || ~isfile( delay_file ) || nFramesDelay ~= nFrames
    origPropDelay = zeros(1, nFrames);
    for idFrame = 1:nFrames
        load(sim_result_files{idFrame},'atmoPath');
        origPropDelay(idFrame) = atmoPath.t(end);
    end
    
    propDelay = origPropDelay;
    if obj.filter_delays && size(sourcePositions,1) > 1
        diffPropDelay = lowPassFilter( diff(origPropDelay), updateRate, 'lowpassiir', updateRate/6, updateRate/3, 1e-4 );
        for idFrame = 1:nFrames
            propDelay(idFrame) = origPropDelay(1) + sum( diffPropDelay(1:idFrame-1) );
        end
    end
    save(delay_file, 'propDelay');
end

%% Determine virtual sources

if force || ~isfile( sim_result_files{end} ) || ~ismember('virtualSource', who('-file', sim_result_files{end}) )
    load(delay_file, 'propDelay');
    virtualSourcePosition = zeros(nFrames, 3);
    
    wb = itaWaitbar(nFrames,'Calculating Virtual Source...');
    for idFrame = 1:nFrames
        load(sim_result_files{idFrame}, 'atmoPath');
        virtualSourcePosition(idFrame,:) = ...
            obj.VirtualSource(atmoPath, receiverPositions(idFrame,:), propDelay(idFrame));
        
        if isvalid(wb); wb.inc(); end
    end
    if isvalid(wb); wb.close(); end
    
    % TODO: Can low pass of VS be improved? Is it necessary to low-pass spreading loss? 
    if size(virtualSourcePosition,1) > 1
        virtualSourcePosition = lowPassFilter(virtualSourcePosition, updateRate, 'lowpassiir', updateRate/6, updateRate/3, 1 );
    end
    
    for idFrame = 1:nFrames
        virtualSource = virtualSourcePosition(idFrame, :);
        smoothedAtmoDelay = propDelay(idFrame);
        save(sim_result_files{idFrame}, 'virtualSource', 'smoothedAtmoDelay', '-append');
    end
end

%% Run urban propagation
if force || ~isfile( sim_result_files{end} ) || ~ismember('urbanPaths', who('-file', sim_result_files{end}) )
    urbanResultFileBackup = obj.urbanPropagation.result_file_path;
    urbanConfFileBackup = obj.urbanPropagation.config_file_path;
    urbanStatsFileBackup = obj.urbanPropagation.runtime_stats_path;
    
    obj.urbanPropagation.result_file_path = fullfile(outputFolder, 'pigeon_urban_paths.json');
    obj.urbanPropagation.config_file_path = fullfile(outputFolder, 'pigeon_project.ini');
    obj.urbanPropagation.runtime_stats_path = fullfile(outputFolder, 'pigeon_project_stats.json');
    
    wb = itaWaitbar(nFrames,'Running urban propagation...');
    for idFrame = 1:nFrames
        load(sim_result_files{idFrame}, 'virtualSource');
        urbanPaths = obj.RunUrbanPropagation( virtualSource, receiverPositions(idFrame, :) );
        save(sim_result_files{idFrame}, 'urbanPaths', '-append');
        if isvalid(wb); wb.inc(); end
    end
    if isvalid(wb); wb.close(); end
    
    obj.urbanPropagation.result_file_path = urbanResultFileBackup;
    obj.urbanPropagation.config_file_path = urbanConfFileBackup;
    obj.urbanPropagation.runtime_stats_path = urbanStatsFileBackup;
end

function [frame_files, delay_file] = GetFileList(outputFolder, nFrames)

frame_files = cell(1, nFrames);
for idFrame = 1:nFrames  
    frame_files{idFrame} = fullfile( outputFolder, ['cau_simulation_frame_',num2str(idFrame),'.mat'] );
end
delay_file = fullfile( outputFolder, 'atmo_prop_delay.mat' );