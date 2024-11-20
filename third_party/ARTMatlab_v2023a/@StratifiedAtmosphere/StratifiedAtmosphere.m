classdef StratifiedAtmosphere
%StratifiedAtmosphere Representing an atmosphere with altitude-dependent weather data
%
% Stores the data to describe parameters of the atmosphere, such as wind,
% temperature and speed of sound. The atmosphere is assumed to be layered
% along the altitude. In which way atmospheric data is dependent on the
% altitude can be adjusted using different profiles for the temperature,
% wind and humidity.
%
% After setting the necessary properties (profiles and constants) a variaty
% of functions can be used to access the speed of sound and the wind
% velocity vector at a specific altitude.

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
    
    
    %% Properties
    properties(Access = public)
        windProfile = 'log';            %Wind profile, see StratifiedAtmosphere.windProfileStrings
        temperatureProfile = 'isa';     %Temperature profile, see StratifiedAtmosphere.temperatureProfileStrings
        humidityProfile = 'constant';   %Humidity profile, see StratifiedAtmosphere.humidityProfileStrings
        
        constTemperature = 293;         %Temperature used in a constant temperature profile
        constStaticPressure = 101325;   %Static pressure used in a constant temperature profile
        
        constRelHumidity = 50;          %Constant Realitive Humidity [%]
        
        constWindDirection = [1 0 0];   %Direction of the wind (a normalized 3D-vector)
        constWindVelocity = 20;         %Wind velocity used in a constant wind profile
        
        frictionVelocity = 0.6;         %Friction velocity used for the logarithmic wind profile
        surfaceRoughness = 0.1;         %Surface roughness used for the logarithmic wind profile

        enableWindPrecalculation = false;       %Switch for precalculating of imported wind data using a vector of 1m resolution
        enableTemperaturePrecalculation = false;%Switch for precalculating of imported temperature data using a vector of 1m resolution
        enableHumidityPrecalculation = false;   %Switch for precalculating of imported humidity data using a vector of 1m resolution
        maxPrecalculationAltitude = 15000;      %Maximum altitude for precalculation of weather data
    end
    properties(Constant = true, Hidden = true)
        ratioSpecHeats = 1.4;           %Ratio of specific heats for air
        airGasConstant = 287.058;       %Gas constant for dry air
        
        karmanConst = 0.4;              %Karman constant
        humidityProfileStrings = {'constant', 'import'};
        temperatureProfileStrings = {'constant', 'isa', 'import'};
        windProfileStrings = {'zero', 'constant', 'log', 'import'};
    end
    
    properties(GetAccess = public, SetAccess = private)
        zMin = 0;       %The lower limit of this atmosphere (default 0)
        %zMax = 12000;   %The upper limit of this atmosphere (default 12000)
    end
    
    properties(GetAccess = public, SetAccess = private, Hidden = true)
        importedZHumidity;
        importedH;
        
        importedZTemperature;
        linFitT;
        importedT;
        importedP0;
        
        importedZWind
        importedVDir;
        importedV;
    end
    
    properties(Access = private, Hidden = true)
        splineT;
        splineP0;
        splineH;
        splinedT;
        splineV;
        splinedV;
%         splineVDir;
%         splinedVDir;
    end
    
    %% Load & Save
    methods
        function store(this, filename)
            % Stores StratifiedAtmosphere into a file using the JSON format
            fID = fopen(filename, 'w');
            try
                fwrite(fID, this.toJSON);
            catch err
                fclose(fID);
                rethrow(err);
            end
            fclose(fID);
        end
    end
    methods(Static = true)
        function atmosphere = load(filename)
            % Loads StratifiedAtmosphere from file using the JSON format
            jsonStr = fileread(filename);
            atmosphere = StratifiedAtmosphere.parseJSON(jsonStr);
        end
        %Parses a json string to a new StratifiedAtmosphere
        atmosphere = parseJSON(jsonStr);
        
        function atmosphere = fromWeatherFile(filename, constWindDirection)
            %Creates a new StratifiedAtmosphere and imports data from a weather
            %file (see import)
            atmosphere = StratifiedAtmosphere();
            if nargin == 2
                atmosphere = atmosphere.import(filename, constWindDirection);
            else
                atmosphere = atmosphere.import(filename);
            end
        end
        function weatherFileToJSON(weatherFileIn, jsonFileOut)
            %Reads a weather file and converts it into a json file
            %representing an StratifiedAtmosphere with imported data.
            atmosphere = StratifiedAtmosphere.load(weatherFileIn);
            atmosphere.store(jsonFileOut);
        end
    end
    
    %% Import Weather Data
    methods
        function this = import(this, filename, constWindDirection)
            %Imports a given weather data file into this atmosphere.
            %Optionally, a constant wind direction can be specified as
            %second parameter (1x3 vector).
            %
            %Imported data is altitude dependent and contains wind vector,
            %temperature, humidity and static pressure along supporting
            %points along the z-axis. Although the latter two are not yet
            %integrated into this class.
            %
            %Later the data between these points is received using spline
            %interpolation. Same holds for the derivatives for wind vector
            %and speed of sound.
            %
            %For importing weather data without this class, see
            %AtmosphericSounding class.
            
            [z, T, vDir, v, humidity, pStat] = AtmosphericSounding.Read(filename);
            if nargin == 3
                assert(isequal(size(constWindDirection), [1 3]), 'Second input must be a 1x3 vector.')
                vDir = repmat(constWindDirection / norm(constWindDirection), size(vDir,1),1);
            end            
            
            this = this.importHumidity(z, humidity);
            this = this.importTemperature(z, T, pStat);
            this = this.importWind(z, v, vDir);
        end
        
        function this = importHumidity(this, z, humidity)
            %Importing a humidity profile using a vector for altitude and
            %a vector with the corresponding relative humidity
            %   Inputs:
            %   z:          Altitude, Nx1 vector [m]
            %   humidity:   Relative humidity, Nx1 vector [%]
            
            if any(z < this.zMin)
                warning('StratifiedAtmosphere: Imported humidity profile has data at negative altitudes.')
            end
            this.humidityProfile = 'import';
            this.importedZHumidity =  z;
            
            this.splineH = spline(z,humidity);
            this.importedH = humidity;
        end
        function this = importTemperature(this, z, T, pStat)
            %Importing a temperature and static pressure profile using a
            %vector for altitude and two vectors with data of corresponding
            %static pressure and temperature
            %   Inputs:
            %   z:      Altitude, Nx1 vector [m]
            %   T:      Temperature [K]
            %   pStat:  Static pressure, Nx1 vector [Pa]
            
            if any(z < this.zMin)
                warning('StratifiedAtmosphere: Imported temperature profile has data at negative altitudes.')
            end
            this.temperatureProfile = 'import';
            this.importedZTemperature =  z;
            
            this.splineP0 = spline(z,pStat);
            this.importedP0 = pStat;
            
            this.linFitT = polyfit(z, T, 1);
            this.splineT = spline(z,T);
            this.splinedT = fnder(this.splineT, 1);
            this.importedT = T;
        end
        function this = importWind(this, z, v, vDir)
            %Importing a wind profile using a arrays for altitude, wind
            %velocity and wind direction vectors.
            %   Inputs:
            %   z: Altitude, Nx1 vector [m]
            %   v: Wind velocity, Nx1 vector [m/s]
            %   vDir: Altitude, Nx3 matrix []
            
            nElem = numel(z);
            if size(vDir, 2) == 2
                vDir = [vDir zeros(size(vDir, 1), 1)];
            end
            if size(vDir, 1) == 1
                vDir = repmat(vDir, nElem, 1);
            end
            
            assert( iscolumn(z), 'z must be a column vector.' )
            assert( iscolumn(v) && numel(v) == nElem, 'v must be a column vector with same size as z.')
            assert( size(vDir, 1) == nElem && size(vDir, 2) == 3, 'vDir must be an Nx3 matrix. N must either be 1 or match the number of altitude values.' )
            assert( all(vDir(:,3) == 0) , 'Wind direction must be purely horizontal.' );
            
            if any(z < this.zMin)
                warning('StratifiedAtmosphere: Imported wind profile has data at negative altitudes.');
            end
            this.windProfile = 'import';
            this.importedZWind =  z;
            
            vTot = bsxfun(@times,v,vDir);
            this.splineV = spline(z, vTot');
            this.splinedV = fnder(this.splineV, 1);
            
            this.importedVDir = vDir;
            this.importedV = v;
        end
    end
    
    %% Set Functions
    methods
        function this = set.humidityProfile(this, profile)
            assert(ischar(profile) && isrow(profile) && any(strcmp(this.humidityProfileStrings, profile)),...
                'Input must be a valid humidity profile string. Type StratifiedAtmosphere.humidityProfileStrings.')
            this.humidityProfile = profile;
        end
        function this = set.windProfile(this, profile)
            assert(ischar(profile) && isrow(profile) && any(strcmp(this.windProfileStrings, profile)),...
                'Input must be a valid wind profile string. Type StratifiedAtmosphere.windProfileStrings.')
            this.windProfile = profile;
        end
        function this = set.temperatureProfile(this, profile)
            assert(ischar(profile) && isrow(profile) && any(strcmp(this.temperatureProfileStrings, profile)),...
                'Input must be a valid temperature profile string. Type StratifiedAtmosphere.temperatureProfileStrings.')
            this.temperatureProfile = profile;
        end
        
        function this = set.constTemperature(this, value)
            assert(isnumeric(value) && numel(value) == 1 && value > 0, 'Input must be a double scalar > 0')
            this.constTemperature = value;
        end
        
        function this = set.constStaticPressure(this, value)
            assert(isnumeric(value) && numel(value) == 1 && value > 0, 'Input must be a double scalar > 0')
            this.constStaticPressure = value;
        end
        
        function this = set.constRelHumidity(this, value)
            assert(isnumeric(value) && numel(value) == 1, 'Input must be a double scalar')
            assert(value >= 0 && value <= 100, 'Input must a percentage between 0 an 100.')
            this.constRelHumidity = value;
        end
        
        
        function this = set.constWindDirection(this, value)
            if iscolumn(value); value = value'; end
            if numel(value) == 2; value = [value 0]; end
            assert(isnumeric(value) && (numel(value) == 3), 'Input must be a double vector with two or three components')
            assert(value(3) == 0, 'Wind direction must be purely horizontal!' )
            if(norm(value) ~=1)
                value = value/norm(value);
                warning('Wind direction must be a normalized vector. Normalizing...');
            end
            this.constWindDirection = value;
        end
        function this = set.constWindVelocity(this, value)
            assert(isnumeric(value) && numel(value) == 1, 'Input must be a double scalar')
            this.constWindVelocity = value;
        end
        function this = set.frictionVelocity(this, value)
            assert(isnumeric(value) && numel(value) == 1, 'Input must be a double scalar')
            this.frictionVelocity = value;
        end
        function this = set.surfaceRoughness(this, value)
            assert(isnumeric(value) && numel(value) == 1 && value > 0, 'Input must be a double scalar > 0')
            this.surfaceRoughness = value;
        end
    end
    
    %% Main Functions - Calculations of Atmospheric Properties
    methods
        function p0 = staticPressure(this, altitude)
            %Calculates the static pressure (Pa) at a given altitude (m).
            %
            %Same as p0()
            switch(this.temperatureProfile)
                case 'constant'
                    p0 = this.constStaticPressure;
                case 'isa'
                    p0 = p0_ISA(this,altitude);
                case 'import'
                    p0 = p0_imp(this, altitude);
                otherwise
                    error('Invalid temperature profile.');
            end
        end
        
        function p = p0(this,altitude)
            %Calculates the static pressure (Pa) at a given altitude (m).
            %
            %Same as staticPressure()
            p = this.staticPressure(altitude);
        end
        
        function alpha = attenuation(this, altitude, f)
            %Calculates the sound attenutation for this atmosphere at a
            %given altitude and frequency bins.
            %
            %Inputs:
            %altitude:  Altitude (m) [1x1 double]
            %f:         Frequency vector (Hz) [1xN / Nx1 double]
            %
            %Outputs:
            %alpha:     Attenuation coefficient (dB/m)
            alpha =  airAttenuationISO(this, altitude, f);
%             alpha =  airAttenuationBASS(this, altitude, f);
        end
        

        function temp = T(this,altitude)
            %Calculates the temperature in Kelvin (K) at the given altitude
            %in meters (m).
            %
            %Same as temperature(altitude)
            switch(this.temperatureProfile)
                case 'constant'
                    temp = this.constTemperature;
                case 'isa'
                    temp = T_ISA(this,altitude);
                case 'import'
                    temp = T_imp(this, altitude);
                otherwise
                    error('Invalid temperature profile.');
            end
        end
        function temp = temperature(this, altitude)
            %Calculates the temperature in Kelvin (K) at the given altitude
            %in meters (m).
            %
            %Same as T(altitude)
            temp = T(this,altitude);
        end
        
        function temp = dT(this,altitude)
            %Calculates the temperature gradient in Kelvin (K/m) at the
            %given altitude in meters (m).
            switch(this.temperatureProfile)
                case 'constant'
                    temp = 0;
                case 'isa'
                    temp = -0.0065;
                case 'import'
                    temp = dT_imp(this, altitude);
                otherwise
                    error('Invalid temperature profile.');
            end
        end
        
        function h = humidity(this, altitude)
            %Returns the relative humidity (percentage between 0 and 100)
            %at a given altitude (m).
            
            switch (this.humidityProfile)
                case 'constant'
                    h = this.constRelHumidity;
                case 'import'
                    h = ppval(this.splineH, altitude);
                otherwise
                    error('Invalid humidity profile.');
            end
        end
        
        function speed = c(this,altitude)
            %Calculates the speed of sound in m/s at the given altitude
            %(m).
            %
            %Same as speedOfSound(altitude)
            speed = sqrt( this.ratioSpecHeats * this.airGasConstant * T(this, altitude) );
        end
        function speed = speedOfSound(this,altitude)
            %Calculates the speed of sound in m/s at the given altitude
            %(m).
            %
            %Same as c(altitude)
            speed = c(this,altitude);
        end
        
        function dSpeed = dc(this,altitude)
            %Calculates the derivative of speed of sound in respect to the
            %altitude at the given altitude.
            %
            %Same as derivedSpeedOfSound(altitude)
            dSpeed = 0.5*sqrt( this.ratioSpecHeats*this.airGasConstant / this.T(altitude) ) * this.dT(altitude); 
        end
        function dSpeed = derivedSpeedOfSound(this,altitude)
            %Calculates the derivative of speed of sound in respect to the
            %altitude at the given altitude.
            %
            %Same as dc(altitude)
            dSpeed = dc(this,altitude);
        end
        
        function wind = v(this,altitude)
            %Calculates the wind velocity vector (m/s) at the given
            %altitude.
            %
            %Same as windVec(altitude)
            switch(this.windProfile)
                case 'zero'
                    wind = [0 0 0];
                case 'constant'
                    wind = this.constWindVelocity * this.constWindDirection;
                case 'log'
                    wind = vLog(this,altitude) * this.constWindDirection;
                case 'import'
                    wind = vImp(this, altitude);
                otherwise
                    error('Invalid wind profile.');
            end
        end
        function wind = windVec(this,altitude)
            %Calculates the wind velocity vector (m/s) at the given
            %altitude.
            %
            %Same as v(altitude)
            wind = v(this,altitude);
        end
        
        function dWind = dv(this,altitude)
            %Calculates the derivative of the wind velocity vector in
            %respect to the altitude at the given altitude.
            %
            %Same as derivedWindVec(altitude)
            switch(this.windProfile)
                case 'zero'
                    dWind = [0 0 0];
                case 'constant'
                    dWind = [0 0 0];
                case 'log'
                    dWind = dvLog(this,altitude) * this.constWindDirection;
                case 'import'
                    dWind = dvImp(this, altitude);
                otherwise
                    error('Invalid wind profile.');
            end            
        end
        function dWind = derivedWindVec(this,altitude)
            %Calculates the derivative of the wind velocity vector in
            %respect to the altitude at the given altitude.
            %
            %Same as dv(altitude)
            dWind = dv(this,altitude);
        end
    end
    
    %% Private Calculation Functions
    methods(Access = private, Hidden = true)
        function temp = T_ISA(~, altitude)
            temp = 288.15;
            if(altitude>0)
                temp = temp - 0.0065*altitude;
            end
        end
        
        function temp = T_imp(this,altitude)
            temp = ppval(this.splineT, altitude);
        end
        
        function p0 = p0_ISA(~,altitude)
            
            p0 = 101325;
            if(altitude>0)
                T0 = 288.15;
                p0 = p0*(1 - 0.0065*altitude/T0)^5.2561;
            end            
        end
        
        function p0 = p0_imp(this,altitude)
            p0 = ppval(this.splineP0, altitude);
        end
        
        function dTemp = dT_imp(this, altitude)
            dTemp = ppval(this.splinedT, altitude);
        end
        
        function wind = vLog(this,altitude)
            if(altitude>this.surfaceRoughness)
                wind = this.frictionVelocity/this.karmanConst * log(altitude / this.surfaceRoughness );
            else
                wind = 0;
            end
        end
        
        function wind = vImp(this,altitude)
            wind = ppval(this.splineV, altitude)';
        end
        
        function dWind = dvLog(this,altitude)
            if(altitude>this.surfaceRoughness)
                dWind =  this.frictionVelocity/this.karmanConst / altitude;
            else
                dWind = 0;
            end
        end
        
        function dWind = dvImp(this,altitude)
            dWind = ppval(this.splinedV, altitude)'; 
        end
        
        A = airAttenuationISO(this, altitude, f);        
    end
end