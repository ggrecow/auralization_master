function atmosphere = parseJSON(jsonStr)
%Parses a JSON formatted string to an StratifiedAtmosphere object

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

%% Init
atmosphere = StratifiedAtmosphere;

%% Decode
jsonStruct = jsondecode(jsonStr);

%% Parse individual profiles
atmosphere = parseHumidity(atmosphere, jsonStruct.humidity_profile);
atmosphere = parseTemperature(atmosphere, jsonStruct.temperature_profile);
atmosphere = parseWind(atmosphere, jsonStruct.wind_profile);


function atmosphere = parseHumidity(atmosphere, jsonStruct)
assert(strcmpi(jsonStruct.type, 'humidity_profile'), ['Trying to parse wrong weather profile type: ' jsonStruct.type ' instead of humidity_profile'])
switch(lower(jsonStruct.class))
    case 'constant'
        atmosphere.humidityProfile = 'constant';
        atmosphere.constRelHumidity = jsonStruct.relative_humidity;
    case 'import'
        atmosphere.humidityProfile = 'import';
        atmosphere = atmosphere.importHumidity(jsonStruct.altitude_values, jsonStruct.humidity_values);

        if isfield(jsonStruct, 'precalculation_enabled')
            atmosphere.enableHumidityPrecalculation = jsonStruct.precalculation_enabled;
        end
        if isfield(jsonStruct, 'precalculation_max_altitude')
            atmosphere.maxPrecalculationAltitude = jsonStruct.precalculation_max_altitude;
        end
    otherwise
        error('Unknown humidity profile class.')
end

function atmosphere = parseTemperature(atmosphere, jsonStruct)
assert(strcmpi(jsonStruct.type, 'temperature_profile'), ['Trying to parse wrong weather profile type: ' jsonStruct.type ' instead of temperature_profile'])
switch(lower(jsonStruct.class))
    case 'constant'
        atmosphere.temperatureProfile = 'constant';
        atmosphere.constStaticPressure = jsonStruct.static_pressure;
        atmosphere.constTemperature = jsonStruct.temperature;
    case 'isa'
        atmosphere.temperatureProfile = 'isa';
    case 'import'
        atmosphere.temperatureProfile = 'import';
        z = jsonStruct.altitude_values;
        pStat = jsonStruct.static_pressure_values;
        T = jsonStruct.temperature_values;
        atmosphere = atmosphere.importTemperature(z, T, pStat);

        if isfield(jsonStruct, 'precalculation_enabled')
            atmosphere.enableTemperaturePrecalculation = jsonStruct.precalculation_enabled;
        end
        if isfield(jsonStruct, 'precalculation_max_altitude')
            atmosphere.maxPrecalculationAltitude = jsonStruct.precalculation_max_altitude;
        end
    otherwise
        error('Unknown temperature profile class.')
end

function atmosphere = parseWind(atmosphere, jsonStruct)
assert(strcmpi(jsonStruct.type, 'wind_profile'), ['Trying to parse wrong weather profile type: ' jsonStruct.type ' instead of wind_profile'])
switch(lower(jsonStruct.class))
    case 'constant'
        atmosphere.windProfile  = 'constant';
        if iscolumn(jsonStruct.wind_vector); jsonStruct.wind_vector = jsonStruct.wind_vector'; end
        atmosphere.constWindVelocity = norm(jsonStruct.wind_vector);
        atmosphere.constWindDirection = jsonStruct.wind_vector / atmosphere.constWindVelocity;
    case 'zero'
        atmosphere.windProfile  = 'zero';
    case 'log'
        atmosphere.windProfile  = 'log';
        atmosphere.frictionVelocity = jsonStruct.friction_velocity;
        atmosphere.surfaceRoughness = jsonStruct.surface_roughness;
        if iscolumn(jsonStruct.wind_direction); jsonStruct.wind_direction = jsonStruct.wind_direction'; end
        atmosphere.constWindDirection = jsonStruct.wind_direction;
    case 'import'
        atmosphere.windProfile  = 'import';
        z = jsonStruct.altitude_values;
        v = jsonStruct.velocity_values;
        vDir = jsonStruct.direction_vectors;
        if size(vDir, 2) == 2
            vDir = [vDir zeros(size(vDir, 1), 1)];
        end
        assert(size(vDir,2) == 3, 'direction_vectors must be a Nx2 (or Nx3) matrix.')
        atmosphere = atmosphere.importWind(z, v, vDir);

        if isfield(jsonStruct, 'precalculation_enabled')
            atmosphere.enableWindPrecalculation = jsonStruct.precalculation_enabled;
        end
        if isfield(jsonStruct, 'precalculation_max_altitude')
            atmosphere.maxPrecalculationAltitude = jsonStruct.precalculation_max_altitude;
        end
    otherwise
        error('Unknown wind profile class.')
end