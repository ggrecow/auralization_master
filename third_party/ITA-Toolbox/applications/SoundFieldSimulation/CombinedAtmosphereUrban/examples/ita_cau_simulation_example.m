%% COMBINED ATMOSPHERIC URBAN SIMULATION
% This interface combines separate simulation tools for atmospheric and
% urban sound propagation. It utilizes a single atmospheric path to model
% the impact of atmospheric refraction and tranlsation on the different
% paths occurring in an urban environment. Based on these paths, the
% overall transfer function is determined.
% 
% The combined simulation requires:
% 1. ARTMatlab (>= v2021b) framework for the atmospheric simulation
% 2. pigeon interface (>= v2021a) for the urban simulation

%% Atmosphere settings
atmos = StratifiedAtmosphere;

atmos.windProfile = 'log';          %String: 'zero', 'constant', 'log'
atmos.temperatureProfile = 'isa';   %String: 'constant', 'isa'
atmos.humidityProfile = 'constant'; %String: 'constant'

atmos.surfaceRoughness = 0.1;       %Surface Roughness for Log Wind Profile [m]
atmos.frictionVelocity = 0.6;       %Friction velocity for Log Wind Profile [m/s]
atmos.constWindDirection = [-1 0 0]; %Normal in wind direction []

atmos.constRelHumidity = 50;        %Constant Realitive Humidity [%]

%% Atmospheric Ray Tracing settings
art = AtmosphericRayTracer;  

art.abortMaxNAdaptations = 50;      %maximum number of ray resolution adaptations
art.abortMinAngleResolution = 1e-5; %minimum angle resolution [°]

art.maxAngleForGeomSpreading = 0.01;%maximimum delta angle of initial direction of neighboring rays used for the calculation of the spreading loss [°]
art.maxPropagationTime = 50;        %maximum propagation time [s]
art.maxReceiverRadius = 0.25;       %accuracy of ray tracer [m]

art.maxReflectionOrder = 0;         %maximum reflection order
art.maxSourceReceiverAngle = 1;     %maximum allowed angle between source and receiver sphere [double between (0, 90) [°]]

%% Urban simulation settings
urbanSim = itaPigeonProject;

urbanSim.MaxDiffractionOrder = 1;   
urbanSim.MaxReflectionOrder = 1;
urbanSim.MaxCombinedOrder = 3;
urbanSim.OnlyNeighbouredEdgeDiffraction = false;
urbanSim.DiffractionOnlyIntoShadowedEdges = false;
urbanSim.FilterNotVisiblePaths = true;
urbanSim.FilterNotVisiblePointToEdge = true;

urbanSim.IntersectionTestResolution = 0.001;
urbanSim.NumIterations = 10;
urbanSim.MaxAccumulatedDiffractionAngle = - 6.2;

urbanSim.LevelDropThreshold = -131.0;
urbanSim.ReflectionPenalty = 0.97;
urbanSim.DiffractionPenalty = 2.5;
urbanSim.result_file_path = [cd,'\results\example.json'];
urbanSim.run_quiet = true;

if ~exist(erase(urbanSim.result_file_path,'example.json'), 'dir')
    mkdir(erase(urbanSim.result_file_path,'example.json'))
end

%% Initialize combined simulation
cauSimulation = itaCauSimulation;

% save settings; not necessary if standard settings defined above are used
cauSimulation.atmosphere = atmos;
cauSimulation.rayTracer = art;
cauSimulation.urbanPropagation = urbanSim;

% chose mode to determine distance between virtual source and receiver
cauSimulation.virtualSourceModus = 'delay'; 
% set path to urban model
cauSimulation.cityModel = which('streetCanyon.skp');

%% Run simulation for source receiver pair
receiverPosition = [17.5,100,1.8];
sourcePosition = [-2000,100,1000];

[atmoPath, virtSource, urbanPaths] = cauSimulation.Run(sourcePosition,receiverPosition);

%% Create TF from result
propagationModel = itaCauPropagation(cauSimulation);
propagationModel.SetSoundPaths(atmoPath, urbanPaths);
tf = propagationModel.TransferFunction();

%% Plot frequency
tf.pf();