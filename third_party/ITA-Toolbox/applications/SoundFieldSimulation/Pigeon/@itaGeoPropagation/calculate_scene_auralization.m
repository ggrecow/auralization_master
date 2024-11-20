function auralization_frame_files = calculate_scene_auralization(obj, pigeon_project, ...
    sim_result_files, auralization_result_folder, freq_vector, force)
%METHOD1 Summary of this method goes here
%   Detailed explanation goes here

arguments
    obj (1,1) itaGeoPropagation
    pigeon_project (1,1) itaPigeonProject
    sim_result_files (1,:) cell
    auralization_result_folder {mustBeTextScalar}
    freq_vector (1,:) double = ita_ANSI_center_frequencies([20 20000], 3)
    force (1,1) logical = false
end

% Init
nFrames = numel(sim_result_files);
obj.sim_prop.orders.combined = pigeon_project.MaxCombinedOrder;
obj.freq_vector = freq_vector.';

% Output files
auralization_frame_files = sim_result_files;

for idx = 1:numel(auralization_frame_files)
    [~,file_name] = fileparts(auralization_frame_files(idx));
    file_name = replace(file_name, 'simulation', 'auralization');
    auralization_frame_files{idx} = fullfile( auralization_result_folder, [file_name,'.mat'] );
end
if ~exist(auralization_result_folder, 'dir')
     mkdir(auralization_result_folder); 
end

% Return if no calculation required
if ~( force || ~isfile( auralization_frame_files{end} ) || ~ismember('paths_update', who('-file', auralization_frame_files{end})) )
    return;
end

% Calculate auralization parameters for all frames
oldUrbanPaths = [];
S = 1; %TODO: Make this a class property
R = 1;
wb = itaWaitbar(nFrames,'Calculating DSP parameters...');
for idFrame = 1:nFrames
    load(sim_result_files{idFrame}, 'urbanPaths');

    % obj.pps_old = oldUrbanPaths;
    % obj.pps = urbanPaths;
    
    disp(['Current Frame ', num2str(idFrame)])
    
    if ~isempty(urbanPaths)

        [spreading_loss, diffraction_gain, delay, air_attenuation, object_interaction, ...
            source_wf_normal, receiver_wf_normal, valid, ~] = ...
            obj.wavefront_coeffs(urbanPaths);
        
        [urbanPaths, hashTable] = ita_propagation_paths_add_identifiers(urbanPaths);
    else
        
    end
    
    paths_update = ita_propagation_va_struct(urbanPaths, oldUrbanPaths, S, R, spreading_loss, diffraction_gain, delay, source_wf_normal, receiver_wf_normal, air_attenuation, object_interaction, [], false);
    
    oldUrbanPaths = urbanPaths;

    save( auralization_frame_files{idFrame}, 'paths_update', 'hashTable' );
    wb.inc();
end
if isvalid(wb); wb.close(); end

end