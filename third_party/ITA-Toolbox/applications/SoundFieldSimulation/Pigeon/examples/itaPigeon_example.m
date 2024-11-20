%% General

pgn = itaPigeonProject();
pgn.geometry_file_path = fullfile( pwd, 'pigeon_test.skp' );
pgn.config_file_path = 'pigeon_test.ini';
pgn.result_file_path = 'pigeon_test.json';
pgn.run_quiet = false;
pgn.export_visualization = true;

emitter_pos = [ 10 10 1.1 ];
sensor_pos = [ 0 0 1.7 ];

paths = pgn.run( emitter_pos, sensor_pos );

fprintf( 'Received %i geometrical propagation paths (see itaGeoPropagation example)\n', numel( paths ) )

%{
pgn.run_gui() % Opens the GUI of pigeon instead
%}
