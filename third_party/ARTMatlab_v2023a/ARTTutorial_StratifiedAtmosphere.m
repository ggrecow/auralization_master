% -------------------------------------------------------------------------
%                        ____  __________  _______
%                       //  / //__   ___/ //  _   |
%                      //  /    //  /    //  /_|  |
%                     //  /    //  /    //  ___   |
%                    //__/    //__/    //__/   |__|
%
% -------------------------------------------------------------------------
%                  ARTMatlab - ITAGeometricalAcoustics
%        (c) Copyright Institute of Technical Acoustics (ITA)
%  This file is part of the ARTMatlab application. Some rights reserved.
% You can find the license in the LICENSE.md file in the ARTMatlab folder.
%--------------------------------------------------------------------------

%% ---ART Tutorial - Stratified Atmosphere---
%  -------------------------------------------
% This is a tutorial on how to setup the StratifiedAtmosphere class using
% different configurations.

%% Initialize class object
atmos = StratifiedAtmosphere;

%% ---Temperature profile---
%  -------------------------
% An altitude-dependent profile for temperature and static pressure.
%
% Note that the static pressure is not used during the Ray Tracing process,
% but for example when calculating the air attenuation based on an
% eigenray.

%% Option 1: Constant (homogeneous) profile
atmos.temperatureProfile = 'constant';
atmos.constTemperature = 293.15; %[K]
atmos.constStaticPressure = 101325; %[Pa]

%% Option 2: International Standard Atmosphere (ISA)
atmos.temperatureProfile = 'isa';

%% Option 3: Import
z = (0:10000)';
T = 293 - 0.65 * z/1000; %linear decreasing temperature
pStat = ones(size(z)) * 101300; %Constant static pressure

%Careful: StratifiedAtmosphere is a value class. When calling a function,
%the result must be stored in your local variable
atmos = atmos.importTemperature(z, T, pStat);
fprintf('When importing the profile type is automatically changed.\natmos.temperatureProfile: ''%s''\n', atmos.temperatureProfile)

%% ---Wind profile---
%  ------------------
% An altitude-dependent profile for wind velocity and horizontal direction.

%% Option 1: Zero
% Medium without any movement.
atmos.windProfile = 'zero';

%% Option 2: Constant
% Constant wind velocity and direction
atmos.windProfile = 'constant';
atmos.constWindVelocity = 4; %[m/s]
atmos.constWindDirection = [1 0];   %Normalized, horizontal vector
atmos.constWindDirection = [1 0 0]; %Passing a zero-valued z-component also works

%% Option 3:Log
% Logarithmic wind profile with constant wind direction
atmos.windProfile = 'log';
atmos.surfaceRoughness = 0.1;       %Surface Roughness for Log Wind Profile [m]
atmos.frictionVelocity = 0.6;       %Friction velocity for Log Wind Profile [m/s]
atmos.constWindDirection = [1 1 0]; %Horizontal vector, will be normalized

%% Option 4: Import
z = (0:10000)';
deltaPhi = 10;
phiWind = ( 0:(deltaPhi/(numel(z)-1)):deltaPhi )';
vDir = [cosd(phiWind) sind(phiWind)]; %Slightly rotating wind direction
v = z * 0.001; %Linearly increasing wind velocity. [m/s]

%Careful: StratifiedAtmosphere is a value class. When calling a function,
%the result must be stored in your local variable
atmos = atmos.importWind(z, v, vDir);
fprintf('When importing, the profile type is automatically changed.\natmos.windProfile: ''%s''\n', atmos.windProfile)


%% ---Humidity profile---
%  ----------------------
% An altitude-dependent profile for relative humidity [%].
% The humidity profile is not used during the Ray Tracing process, but for
% example when calculating the air attenuation based on an eigenray.

%% Option 1: Constant
atmos.humidityProfile = 'constant';
atmos.constRelHumidity = 50; % [%]

%% Option 2: Import
z = (0:10000)';
h = 50 + z * 0.001; %Linearly increasing relative humidity. [%]

%Careful: StratifiedAtmosphere is a value class. When calling a function,
%the result must be stored in your local variable
atmos = atmos.importHumidity(z, h);
fprintf('When importing, the profile type is automatically changed.\natmos.humidityProfile: ''%s''\n', atmos.humidityProfile)



%% ---Import atmospheric soundings---
%  ----------------------------------
% The ART framework is capable to import weather data from atmospheric
% soundings provided on:
% http://weather.uwyo.edu/upperair/sounding.html
% 
% Instructions:
% On that page, select a timestamp, click on a weather station and a new
% tab containing the weather data will open. Select the full text (Ctrl + A),
% copy it (Ctrl + C) and save the content in a .txt-file. Now you can
% import from that file as shown below.
% 
% Important note:
% During the import, the altitude data is shifted, so that z=0 refers to
% the first measurement at the elevation of weather station.

%% Define your file here
soundingFilepath = fullfile('@AtmosphericSounding', 'AtmosphericSoundingExample.txt');

%% Create object directly from file
atmos = StratifiedAtmosphere.fromWeatherFile( soundingFilepath );
fprintf( 'Temperature at 100m: %d K\n', atmos.temperature(100) );

%% Read data into existing object
atmos = atmos.import( soundingFilepath );
fprintf( 'Wind velocity at 100m: %d m/s\n', norm(atmos.windVec(100)) );

%% Apply constant wind direction 
% Optionally, the wind direction can be overwritten with a user-defined
% constant direction
atmos = StratifiedAtmosphere.fromWeatherFile( soundingFilepath, [0 -1 0] );
atmos = atmos.import( soundingFilepath, [0 -1 0] );
vDir = atmos.windVec(100) / norm(atmos.windVec(100));
fprintf( 'Wind direction at 100m: [%d %d %d]\n', vDir(1), vDir(2), vDir(3) );



%% ---Acceleration for imported profiles---
%  ----------------------------------------
% When working with imported weather data, the data is interpolated during
% the ray tracing process using a cubic spline approach. This might add a
% certain computational load. In order to speed-up the process, the data
% can be precomputed using an altitude-grid with 1m resolution. Then, the
% data will be acessed using a nearest-neighbor approach.
% NOTE:
% This is only done during the actual ray tracing process (C++ code), not
% when evaluating the weather parameters from Matlab.

% Enable the precalculation for temperature, wind and humidity, respectively.
atmos.enableTemperaturePrecalculation = true;
atmos.enableWindPrecalculation = true;
atmos.enableHumidityPrecalculation = true;

% The maximum altitude considered for the precalculation.
atmos.maxPrecalculationAltitude = 15000;

% Important note:
% The profile will be assumed constant above this maximum altitude, so it
% is advised to select an altitude that covers all potential sound paths.
