function jsonStr = toJSON(atmosphere)
%Converts StratifiedAtmosphere class object to JSON formatted string

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

humidityStruct = getHumidityJSONStruct(atmosphere);
temperatureStruct = getTemperatureJSONStruct(atmosphere);
windStruct = getWindJSONStruct(atmosphere);

jsonStruct.humidity_profile = humidityStruct;
jsonStruct.temperature_profile = temperatureStruct;
jsonStruct.wind_profile = windStruct;

jsonStr = jsonencode(jsonStruct);

function out = getHumidityJSONStruct(this)
out.type = 'humidity_profile';
switch this.humidityProfile
    case 'constant'
        out.class = 'constant';
        out.relative_humidity = this.constRelHumidity;
    case 'import'
        out.class = 'import';
        out.altitude_values = this.importedZHumidity;
        out.humidity_values = this.importedH;
        out.precalculation_enabled = this.enableHumidityPrecalculation;
        out.precalculation_max_altitude = this.maxPrecalculationAltitude;
    otherwise
        error('Humidity profile %s is not supported for JSON export.', char(this.humidityProfile))
end

function out = getTemperatureJSONStruct(this)
out.type = 'temperature_profile';
switch this.temperatureProfile
    case 'constant'
        out.class = 'constant';
        out.static_pressure = this.constStaticPressure;
        out.temperature = this.constTemperature;
    case 'isa'
        out.class = 'isa';
    case 'import'
        out.class = 'import';
        out.altitude_values = this.importedZTemperature;
        out.static_pressure_values = this.importedP0;
        out.temperature_values = this.importedT;
        out.precalculation_enabled = this.enableTemperaturePrecalculation;
        out.precalculation_max_altitude = this.maxPrecalculationAltitude;
    otherwise
        error('Temperature profile %s is not supported for JSON export.', char(this.tempProfile))
end

function out = getWindJSONStruct(this)
out.type = 'wind_profile';
switch this.windProfile
    case 'zero'
        out.class = 'zero';
    case 'constant'
        out.class = 'constant';
        out.wind_vector = this.windVec(0);
    case 'log'
        out.class = 'log';
        out.friction_velocity = this.frictionVelocity;
        out.surface_roughness = this.surfaceRoughness;
        out.wind_direction = this.constWindDirection;
    case 'import'
        out.class = 'import';
        out.altitude_values = this.importedZWind;
        out.velocity_values = this.importedV;
        out.direction_vectors = this.importedVDir(:, 1:2);
        out.precalculation_enabled = this.enableWindPrecalculation;
        out.precalculation_max_altitude = this.maxPrecalculationAltitude;
    otherwise
        error('Wind profile %s is not supported for JSON export.', char(this.windProfile))
end