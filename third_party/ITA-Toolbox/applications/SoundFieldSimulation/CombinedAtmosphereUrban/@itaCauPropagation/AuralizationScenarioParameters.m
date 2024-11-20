function [auralization_frame_files] = AuralizationScenarioParameters(obj, sim_result_files, outputFolder, freqVector, force)
%AURALIZATIONSCENARIOPARAMETERS Calculates parameters used for auralization
% For all results, the second dimension refers to the index of the
% respective sound paths. The first dimension, refers to frequency for
% frequency-dependent parameters.

if nargin < 4; freqVector = ita_ANSI_center_frequencies([20 20000], 3); end
if nargin < 5; force = false; end

%Ensure column vector for frequencies
if isrow(freqVector); freqVector = freqVector'; end

%% Init
nFrames = numel(sim_result_files);

%% Output files
auralization_frame_files = GetFileList(outputFolder, nFrames);
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end

%% Return if no calculation required
if ~( force || ~isfile( auralization_frame_files{end} ) || ~ismember('paths_update', who('-file', auralization_frame_files{end})) )
    return;
end

%% Calculate auralization parameters for all frames
oldUrbanPaths = [];
S = 1; %TODO: Make this a class property
R = 1;
wb = itaWaitbar(nFrames,'Calculating DSP parameters...');
for idFrame = 1:nFrames
    load(sim_result_files{idFrame}, 'atmoPath', 'urbanPaths', 'smoothedAtmoDelay');
    
    obj.SetSoundPaths(atmoPath, urbanPaths);
    [airAbsorption, objectInteraction, spreadingLoss, diffractionGain, delay, sourceWFNormals, receiverWFNormals]...
        = obj.AuralizationParameters(freqVector, smoothedAtmoDelay);


    [urbanPaths, hashTable] = ita_propagation_paths_add_identifiers(urbanPaths);
    paths_update = ita_propagation_va_struct(urbanPaths, oldUrbanPaths, S, R, spreadingLoss, diffractionGain, delay, sourceWFNormals, receiverWFNormals, airAbsorption, objectInteraction, [], false);
    oldUrbanPaths = urbanPaths;
    
    save( auralization_frame_files{idFrame}, 'paths_update', 'hashTable' );
    wb.inc();
end
if isvalid(wb); wb.close(); end


function [auralization_frame_files] = GetFileList(outputFolder, nFrames)

auralization_frame_files = cell(1, nFrames);
for idFrame = 1:nFrames  
    auralization_frame_files{idFrame} = fullfile( outputFolder, ['cau_auralization_frame_',num2str(idFrame),'.mat'] );
end
