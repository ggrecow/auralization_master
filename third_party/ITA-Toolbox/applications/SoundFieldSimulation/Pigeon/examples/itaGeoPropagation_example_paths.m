%% General

% initialize GeoPropSim
geoPropSim = itaGeoPropagation();

% load directivity
% directivity_id = geoPropSim.load_directivity( 'Genelec8020.ir.daff', 'Genelec8020' );

%% run simulation


% Single path
geoPropSim.load_paths( 'ppa_example_paths.json' );

geoPropSim.pps = geoPropSim.pps( 12 ); % Pick a path
pps1TF = itaAudio();
pps1TF.freqData = geoPropSim.run();


% Multiple paths
geoPropSim.load_paths( 'ppa_example_paths_2.json' );

geoPropSim.pps = geoPropSim.pps( [ 1:83 85:181 182 183 185:249 ] ); % Skip paths that lead to errors

pps2TF = itaAudio();
pps2TF.freqData = geoPropSim.run();
