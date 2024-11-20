classdef AtmosphericSounding
%AtmosphericSounding Static class for dealing with weather data (sounding) files
%   Files originate from:
%   University of Wyoming, College of Engineering, Department of
%   Atmospheric Science. Soundings.
%   http://weather.uwyo.edu/upperair/sounding.html

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
    
    methods
        function obj = AtmosphericSounding()
            error('Static class. Constructor is diabled.')
        end
    end
    methods(Static = true)
        function [z, T, vDir, v, humidity, pStat] = Read(filename)
            %Imports altitude-dependent weather data from an atmospheric sounding file.
            %
            %Input:
            %filename:  Name of data file
            %
            %Outputs:
            %z:         Supporting points of altitude [m]
            %T:         Temperature [K]
            %vDir:      Wind direction
            %v:         Wind velocity [m/s]
            %humidity:  Relative humidity [%]
            %pStat:     Static pressure [Pa]
            
            %% Extract Relevant Text-Lines
            % filename = 'WeatherData/FAUPUptington00z04Nov2016.txt';
            
            txt = readTextLinewise(filename);
            [dataTxt, stationInfoTxt] = seperateDataFromStationInfo(txt);
            
            dataNum = dataTxtToNum(dataTxt);
            statElevation = str2double( getStationProperty(stationInfoTxt, 'Station elevation') );
            
            %% Extract Single Data Vectors
            z = dataNum(:, 2) - statElevation;
            T = dataNum(:, 3) + 273.15; %[°C] -> [K]
            vAzimuth = dataNum(:, 7);
            %Conversion from wind degree (0 = wind from north, 90 = wind from east)
            %to wind vector using East-North-Up system (x = wind towards east, y = wind towards north).
            vDir = [-sind(vAzimuth) -cosd(vAzimuth) zeros(size(vAzimuth))];
            v = dataNum(:, 8)*0.514444444; %[knots] -> [m/s]
            
            humidity = dataNum(:, 5);
            pStat = dataNum(:, 1) * 100; %[hPa] -> [Pa]
        end
    end
end

