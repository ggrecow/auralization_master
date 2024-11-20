function atmos = initializeAtmosphere(windDirection)

atmos = StratifiedAtmosphere;
atmos.windProfile = 'log';          %String: 'zero', 'constant', 'log'
atmos.temperatureProfile = 'isa';   %String: 'constant', 'isa'
atmos.humidityProfile = 'constant'; %String: 'constant'

atmos.surfaceRoughness = 0.1;       %Surface Roughness for Log Wind Profile [m]
atmos.frictionVelocity = 0.6;       %Friction velocity for Log Wind Profile [m/s]
atmos.constWindDirection = [-1 0 0]; %Normal in wind direction []

atmos.constRelHumidity = 50;        %Constant Realitive Humidity [%]

if nargin
    if isnumeric(windDirection) && isscalar(windDirection)
        windDirection = WindAzimuthToDirection(windDirection);
    end
    assert(isnumeric(windDirection) && numel(windDirection) == 3, 'Input must be a 1x3 normal vector or an azimuth angle.')

    atmos.constWindDirection = windDirection;
end