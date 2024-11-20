classdef AtmosphericRayTracer < handle
%AtmosphericRayTracer Class to interface the AtmosphericRayTracing framework (part of ITAGeometricalAcoustics)

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

    properties(Access = public)
        %% Ray Tracing Settings
        %--General--
        
        maxPropagationTime = 15;                    %Maximum propagation time of rays [s]
        integrationTimeStep = 0.1;                  %Base time-step for integration of rays [s]
        solverMethod = 'runge-kutta';               %Method for solving the ODE ('runge-kutta' or 'euler')
        multiThreadingActive = true;                %If enabled, tracing multiple rays is done using multi-threading (openMP)
        
        %--Adaptive integration-- (EXPERT SETTINGS)
        
        adaptiveIntegrationActive = true;           %Activate adaptive integration of ray, adapting time step (boolean)
        adaptiveIntegrationMaxError = 0.015;        %Threshold for decreasing time step []
        adaptiveIntegrationUncritivalError = 0.005; %Threshold for re-increasing time step []
        adaptiveIntegrationMaxLevel = 31;           %Maximum adaptation depth for adaptive integration [] (integer)
        
        
        %% Eigenray Search Settings
        %--General--
        
        %maxPropagationTime;                %Maximum propagation time of rays [s] (see above)
        
        maxReflectionOrder = 1;             %Maximum number of considered reflections (=number of eigenrays) [] (integer)
        abortOnReceiverDistIncrease = true; %If enabled, ray tracing will be aborted as soon as ray receiver distance increases (boolean)
        
        %--Adaptive ray zooming--
        
        maxReceiverRadius = 1;              %Maximum value for receiver radius [m]
        maxSourceReceiverAngle = 1;         %Maximum allowed angle between source and receiver sphere [double between (0, 90) [°]]
        maxAngleForGeomSpreading = 0.01;    %Maximimum delta angle of initial direction of neighboring rays used for the calculation of the spreading loss [°]
        
        abortMaxNAdaptations = 30;          %Abort criterion for adaptive ray zooming: maximum number of iterations [] (integer)
        abortMinAngleResolution = 0.001;    %Abort criterion for adaptive ray zooming: minimum angular resolution of rays [°]
        
        
        %--Advanced ray zooming-- (EXPERT SETTINGS)
        
        advancedRayZoomingActive = false;   %If activated this leads to additional run-time improvement but may also cause the eigenray search to fail
        advancedRayZoomingThreshold = 0.1;  %Threshold for advanced ray zooming [double between [0, 2], 0 = always use advanced method, 2 = off]
        
    end
    
    %% Static Conversion Functions
    methods(Static = true)
        function n = AnglesToNormalVectors(thetaDeg, phiDeg, bLaunchAngleConvention)
            %Converts the angles theta (elevation) and phi (azimuth) to an Nx3 matrix representing the corresponding direction
            %   Inputs (defaults):
            %   thetaDeg:   Either vector or matrix of elevation angles [°]
            %   phiDeg:     Either vector or matrix of azimuth angles [°]
            %   bLaunchAngleConvention (true):  Switch for convention of elevation angle
            %   
            %   Note: ART uses the convention of spherical coordinates for
            %   the elevation (north/south pole 0°/180°). In literature,
            %   the initial elevation is often referred to as launch angle
            %   which has a different convention (north/south pole 90°/-90°).
            
            if nargin < 3; bLaunchAngleConvention = false; end
            if bLaunchAngleConvention
                thetaDeg = 90 - thetaDeg;
            end
            
            sinTheta = sind(thetaDeg(:));
            nx = sinTheta.*cosd(phiDeg(:));
            ny = sinTheta.*sind(phiDeg(:));
            nz = cosd(thetaDeg(:));
            n = [nx, ny, nz];
        end
    end
    
    %% Public Interface Methods
    methods(Access = public)
        function obj = AtmosphericRayTracer()
        end
        
        function [eigenrays, totalNumRaysTraced] = FindEigenrays(obj, atmosphere, source, receiver)
            %Finds eigenrays giving atmosphere, source and receiver position
            % As second optional return value, this function provides the
            % total number of traced rays for each eigenray
            if isa(source, 'itaCoordinates'); source = source.cart; end
            if isa(receiver, 'itaCoordinates'); receiver = receiver.cart; end
            
            assert(isa(atmosphere, 'StratifiedAtmosphere'), 'First input must be an object of class StratifiedAtmosphere')
            assert(isnumeric(source) && numel(source) == 3,...
                'Second input must either be a 3-element numeric vector or an itaCoordinates object with one point.')
            assert(isnumeric(receiver) && numel(receiver) ==3,...
                'Third input must either be a 3-element numeric vector or an itaCoordinates object with one point.')
            assert(source(3) >= 0, 'Source positition must be above the ground (z >=0 )')
            assert(receiver(3) >= 0, 'Receiver positition must be above the ground (z >=0 )')
            
            jsonStrAtmosphere = atmosphere.toJSON();
            [jsonStrEigenrays, totalNumRaysTraced] = ARTMatlab('FindEigenrays', jsonStrAtmosphere, source, receiver, obj.GetEigenraySettings(), obj.GetRayTracingSettings());
            eigenrays = AtmosphericRay.parseJSON(jsonStrEigenrays);
        end
        
        function rays = TraceRays(obj, atmosphere, source, rayNormalVectors)
            %Trace rays giving atmosphere, source position and initial ray directions as normal vectors
            if isa(source, 'itaCoordinates'); source = source.cart; end
            if isa(rayNormalVectors, 'itaCoordinates'); rayNormalVectors = rayNormalVectors.cart; end
            
            assert(isa(atmosphere, 'StratifiedAtmosphere'), 'First input must be an object of class StratifiedAtmosphere')
            assert(isnumeric(source) && numel(source) == 3,...
                'Second input must either be a 3-element numeric vector or an itaCoordinates object with one point.')
            assert(isnumeric(rayNormalVectors) && size(rayNormalVectors,2) == 3 && numel(rayNormalVectors) > 0,...
                'Third input must either be a Nx3 numeric vector or an itaCoordinates object with at least one point.')
            assert(source(3) >= 0, 'Source positition must be above the ground (z >=0 )')
            
            
            jsonStrAtmosphere = atmosphere.toJSON();
            jsonStrEigenrays = ARTMatlab('TraceRays', jsonStrAtmosphere, source, rayNormalVectors, obj.GetRayTracingSettings());
            rays = AtmosphericRay.parseJSON(jsonStrEigenrays);
        end
        
    end
    
    methods(Static = true)
        function Info()
            %Displays information about ARTMatlab including its version
            
            %Adjust blanks depending version string length
            nCharsPerRow = 59;
            versionStr = ARTMatlab('Version');
            nMainChars = 12;
            nBlanksAtEnd = nCharsPerRow-nMainChars-numel(versionStr);
            
            disp('* ------------------------------------------------------- *');
            disp('*              _______   _______    __________            *')
            disp('*             //  _   | ||  _   |  //__   ___/            *');
            disp('*            //  /_|  | || |_|  |    //  /                *');
            disp('*           //  ___   | ||  __  \   //  /                 *');
            disp('*          //__/   |__| ||_|  \__\ //__/                  *');
            disp('*                      ARTMatlab                          *');
            disp('* ------------------------------------------------------- *');
            disp('* Interface to the Atmospheric Ray Tracing framework      *');
            disp('* Part of ITAGeometricalAcoustics                         *');
            disp('* https://git.rwth-aachen.de/ita/ITAGeometricalAcoustics/ *')
            disp(['* Version: ' versionStr blanks(nBlanksAtEnd) '*']);
            disp('* ------------------------------------------------------- *');
        end
    end
    
    %% Hidden Interface Methods
    methods(Access = public, Hidden = true)
        function measuredRunTimes = BenchmarkEigenraySearch(obj, atmosphere, source, receiver, numRuns)
            if isa(source, 'itaCoordinates'); source = source.cart; end
            if isa(receiver, 'itaCoordinates'); receiver = receiver.cart; end
            
            assert(isa(atmosphere, 'StratifiedAtmosphere'), 'First input must be an object of class StratifiedAtmosphere')
            assert(isnumeric(source) && numel(source) == 3,...
                'Second input must either be a 3-element numeric vector or an itaCoordinates object with one point.')
            assert(isnumeric(receiver) && numel(receiver) ==3,...
                'Third input must either be a 3-element numeric vector or an itaCoordinates object with one point.')
            assert(isnumeric(numRuns) && isscalar(numRuns), 'Fourth input must be a numeric scalar representing the number of simulation runs.')
            
            jsonStrAtmosphere = atmosphere.toJSON();
            measuredRunTimes = ARTMatlab('BenchmarkEigenraySearch', jsonStrAtmosphere, source, receiver, numRuns, obj.GetEigenraySettings(), obj.GetRayTracingSettings());
        end
    end
    
    %% Parsing settings
    methods(Access = private)
        function settings = GetRayTracingSettings(obj)
            settings.maxPropagationTime = obj.maxPropagationTime;
            settings.integrationTimeStep = obj.integrationTimeStep;
            settings.solverMethod = obj.solverMethod;
            settings.multiThreadingActive = obj.multiThreadingActive;
            settings.adaptiveIntegrationActive = obj.adaptiveIntegrationActive;
            settings.adaptiveIntegrationMaxError = obj.adaptiveIntegrationMaxError;
            settings.adaptiveIntegrationUncritivalError = obj.adaptiveIntegrationUncritivalError;
            settings.adaptiveIntegrationMaxLevel = obj.adaptiveIntegrationMaxLevel;
        end
        
        function settings = GetEigenraySettings(obj)
            settings.maxPropagationTime = obj.maxPropagationTime;
            settings.maxReflectionOrder = obj.maxReflectionOrder;
            settings.abortOnReceiverDistIncrease = obj.abortOnReceiverDistIncrease;
            settings.maxReceiverRadius = obj.maxReceiverRadius;
            settings.maxSourceReceiverAngle = obj.maxSourceReceiverAngle;
            settings.maxAngleForGeomSpreading = obj.maxAngleForGeomSpreading;
            settings.advancedRayZoomingActive = obj.advancedRayZoomingActive;
            settings.advancedRayZoomingThreshold = obj.advancedRayZoomingThreshold;
            settings.abortMaxNAdaptations = obj.abortMaxNAdaptations;
            settings.abortMinAngleResolution = obj.abortMinAngleResolution;
        end
    end
end