function generate_config( obj, source_pos, receiver_pos )
% GENERATE_CONFIG This function generates the .ini file of the city path algorithm

[ output_folder, file_base_name ] = fileparts( obj.config_file_path );

if isempty( output_folder )
    output_folder = pwd;
else
    if ~exist( output_folder, 'dir' )
        mkdir( output_folder )
    end
end

visualizationPaths = [];

for idx = 1:length(obj.visualizationPaths)
       
%    visualizationPaths{idx} = varargin{idx};  
   visualizationPaths{idx} = strrep(obj.visualizationPaths{idx},'\','\\');
    
end

receiverString = ['SensorPos = ',num2str(receiver_pos(1)),', ',num2str(receiver_pos(2)),', ',num2str(receiver_pos(3)),'\n' ];
sourceString =  ['EmitterPos = ',num2str(source_pos(1)),', ',num2str(source_pos(2)),', ',num2str(source_pos(3)),'\n' ];


%% save old content of file as backup

if obj.saveBackupConfig
    
    config_path_backup = fopen( fullfile( output_folder, strcat( file_base_name, '_backup.ini' ) ), 'w+' );

    try
        txt = fileread([obj.outFilePath,obj.resultName,'.ini']);
    catch
        txt = ''; %.ini File not created yet
    end

    % change '\' to '\\'
    txt = strrep(txt,'\','\\');
    fprintf( config_path_backup, txt );
    fclose( config_path_backup );
    
end


%% create new .ini file
file = fopen( obj.config_file_path, 'w' );

% go through file one by one and write the necessary output
fprintf(file,'[pigeon:scene]\n');
fprintf(file,' \n');

fprintf(file, 'GeometryFilePath = %s\n', obj.geometry_file_path );
fprintf(file, 'OutputFilePath = %s\n', obj.result_file_path );
fprintf(file,receiverString);
fprintf(file,sourceString);
fprintf(file,' \n');

fprintf(file,'[pigeon:config] \n');
fprintf(file,' \n');

fprintf(file,['MaxDiffractionOrder = ',num2str(obj.MaxDiffractionOrder),' \n']);
fprintf(file,['MaxReflectionOrder = ',num2str(obj.MaxReflectionOrder),' \n']);
fprintf(file,['MaxCombinedOrder = ',num2str(obj.MaxCombinedOrder),' \n']);
fprintf(file,' \n');

fprintf(file,['OnlyNeighbouredEdgeDiffraction = ',obj.OnlyNeighbouredEdgeDiffraction,' \n']);
fprintf(file,['DiffractionOnlyIntoShadowedEdges = ',obj.DiffractionOnlyIntoShadowedEdges,' \n']);
fprintf(file,['FilterNotVisiblePaths = ',obj.FilterNotVisiblePaths,' \n']);
fprintf(file,['FilterEmitterToEdgeIntersectedPaths = ',obj.FilterNotVisiblePointToEdge,' \n']);
fprintf(file,['FilterSensorToEdgeIntersectedPaths = ',obj.FilterNotVisiblePointToEdge,' \n']);
fprintf(file,['IntersectionTestResolution = ',num2str(obj.IntersectionTestResolution),' \n']);
fprintf(file,['NumIterations = ',num2str(obj.NumIterations),' \n']);
fprintf(file,' \n');

fprintf(file,['MaxAccumulatedDiffractionAngle = ',num2str(obj.MaxAccumulatedDiffractionAngle),' \n']);
fprintf(file,' \n');

fprintf(file,['LevelDropThreshhold = ',num2str(obj.LevelDropThreshold),' \n']);
fprintf(file,['ReflectionPenalty = ',num2str(obj.ReflectionPenalty),' \n']);
fprintf(file,['DiffractionPenalty = ',num2str(obj.DiffractionPenalty),' \n']);
fprintf(file,' \n');

fprintf(file,['ExportVisualisation = ',obj.export_visualization,' \n']);
fprintf(file,' \n');

% fprintf(file,'[pigeon:visualization:pigeonVizLayer1] \n');

for idx=1:length(visualizationPaths)
    
    fprintf(file,['[pigeon:visualization:pigeonVizLayer',num2str(idx),'] \n']);    
    fprintf(file,' \n');
    fprintf(file,['path',num2str(1),' = ',visualizationPaths{idx},' \n']);   
    fprintf(file,' \n');
        
end

fclose(file);

end

