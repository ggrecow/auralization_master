function export_config( obj, s, r )

[ output_folder, file_base_name ] = fileparts( obj.config_file_path );

ini = IniConfig();
ini.AddSections( { 'pigeon:scene', 'pigeon:config' } );

ini.AddKeys( 'pigeon:scene', 'GeometryFilePath', obj.geometry_file_path );
ini.AddKeys( 'pigeon:scene', 'EmitterPos', s );
ini.AddKeys( 'pigeon:scene', 'SensorPos', r );
ini.AddKeys( 'pigeon:scene', 'OutputFilePath', obj.result_file_path );

ini.AddKeys( 'pigeon:config', 'ExportRuntimeStatistics', double( obj.export_runtime_statistics ) );
ini.AddKeys( 'pigeon:config', 'RuntimeStatisticsFilePath', obj.runtime_stats_path );

ini.AddKeys( 'pigeon:config', 'MaxDiffractionOrder', obj.MaxDiffractionOrder );
ini.AddKeys( 'pigeon:config', 'MaxReflectionOrder', obj.MaxReflectionOrder );
ini.AddKeys( 'pigeon:config', 'MaxCombinedOrder', obj.MaxCombinedOrder );

ini.AddKeys( 'pigeon:config', 'OnlyNeighbouredEdgeDiffraction', double( obj.OnlyNeighbouredEdgeDiffraction ) );
ini.AddKeys( 'pigeon:config', 'DiffractionOnlyIntoShadowedEdges', double( obj.OnlyNeighbouredEdgeDiffraction ) );
ini.AddKeys( 'pigeon:config', 'FilterNotVisiblePathsBetweenEdges', double( obj.FilterNotVisiblePointToEdge ) );
ini.AddKeys( 'pigeon:config', 'FilterEmitterToEdgeIntersectedPaths', double( obj.FilterEmitterToEdgeIntersectedPaths ) );
ini.AddKeys( 'pigeon:config', 'FilterSensorToEdgeIntersectedPaths', double( obj.FilterSensorToEdgeIntersectedPaths ) );
ini.AddKeys( 'pigeon:config', 'FilterNotVisiblePaths', double( obj.FilterNotVisiblePaths ) );
ini.AddKeys( 'pigeon:config', 'IntersectionTestResolution', obj.IntersectionTestResolution );
ini.AddKeys( 'pigeon:config', 'NumIterations', obj.NumIterations );

ini.AddKeys( 'pigeon:config', 'MaxAccumulatedDiffractionAngle', obj.MaxAccumulatedDiffractionAngle );
ini.AddKeys( 'pigeon:config', 'LevelDropThreshold', obj.LevelDropThreshold );
ini.AddKeys( 'pigeon:config', 'ReflectionPenalty', obj.ReflectionPenalty );
ini.AddKeys( 'pigeon:config', 'DiffractionPenalty', obj.DiffractionPenalty );

ini.AddKeys( 'pigeon:config', 'ExportVisualisation', double( obj.export_visualization ) );

for idx = 1:numel( obj.visualizationPaths )
    secname = sprintf( 'pigeon:visualization:pigeonVizLayer%00i', idx );
    ini.AddSections( secname );
    ini.AddKeys( secname, sprintf( 'path%00i', idx ), obj.visualizationPaths{ idx } );
end


%% Export INI file

if isempty( output_folder )
    output_folder = pwd;
else
    if ~exist( output_folder, 'dir' )
        mkdir( output_folder );
    end
end

export_file_path = fullfile( output_folder, strcat( file_base_name, '.ini' ) ) ;

if obj.saveBackupConfig && exist( export_file_path, 'file' ) ~= 0
    export_file_path_backup = fullfile( output_folder, stract( file_base_name, '_backup.ini' ) ) ;
    copyfile( export_file_path, export_file_path_backup );
end

ini.WriteFile( export_file_path );


end
