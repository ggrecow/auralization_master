classdef AtmosphericPropagation
%AtmosphericPropagation Allows to derive auralization parameters from eigenrays using suitable propagation models
% It also is able to generate a transfer function based on a set of eigenrays.

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
        atmosphere = StratifiedAtmosphere; %StratifiedAtmosphere object
        frequencyVector = ita_ANSI_center_frequencies([20 20000], 3)'; %Frequency vector for the calculation of the transfer function
        groundReflectionFactor = 1; %Reflection factor used for ground reflection
    end
    
    %% Set
    methods
        function obj = set.atmosphere(obj, atmos)
            assert(isa(atmos, 'StratifiedAtmosphere'),...
                [mfilename('class') ': Input for atmosphere must be an StratifiedAtmosphere object']);
            
            obj.atmosphere = atmos;
        end
        function obj = set.frequencyVector(obj, f)
            assert(isnumeric(f) && isvector(f),...
                [mfilename('class') ': Input for frequencyVector must be numeric vector']);
            
            if isrow(f); f = f'; end
            obj.frequencyVector = f;
        end
        function obj = set.groundReflectionFactor(obj, R)
            assert(isnumeric(R) && ( isscalar(R) || isvector(R) && numel(R) == numel(obj.frequencyVector) ),...
                [mfilename('class') ': Input for groundReflectionFactor must be scalar or have same size as frequencyVector']); %#ok<MCSUP>
            
            if isrow(R); R = R.'; end
            obj.groundReflectionFactor = R;
        end
    end
    
    %% Auralization parameters    
    methods(Static = true)
        function spreadingLoss = SpreadingLoss(eigenrays)
            %Returns the spreading loss of given eigenrays as array of size 1 x nRays.
            spreadingLoss = ones(1, numel(eigenrays));
            for idx = 1:numel(eigenrays)
                if eigenrays(idx).numPoints == 0; continue; end
                spreadingLoss(idx) = eigenrays(idx).spreadingLoss;
            end
        end
        function delay = PropagationDelay(eigenrays)
            %Returns the propagation delay of given eigenrays as array of size 1 x nRays.
            delay = zeros(1, numel(eigenrays));
            for idx = 1:numel(eigenrays)
                if eigenrays(idx).numPoints == 0; continue; end
                delay(idx) = eigenrays(idx).t(end);
            end
        end
    end
    
    methods (Access = public)
        function alpha = AirAttenuation(obj, eigenrays)
            %Calculates the air attenuation factor for given eigenrays
            %
            %Calculation is done by integrating over the attenuation
            %coefficient defined in ISO 9613-1 [dB/m] along the ray path
            %numerically for each frequency. The output has a size of
            %nFreqs x nRays.
            %
            %Outputs:
            %alpha:     Air attenuation factor [0...1] (dimensionless)
            %           Hint: This is no dB value but an actual factor
            
            alphaDB = zeros(numel(obj.frequencyVector), numel(eigenrays));
            
            for idxRay = 1:numel(eigenrays)
                eigenray = eigenrays(idxRay);
                if eigenray.numPoints == 0; continue; end
                
                deltaR = vecnorm(diff(eigenray.r.cart), 2, 2);
                for idxPoint = 1:(eigenray.numPoints-1)
                    z = abs( eigenray.r.z(idxPoint) );
                    alphaDB(:, idxRay) = alphaDB(:, idxRay) + obj.atmosphere.attenuation(z, obj.frequencyVector) * deltaR(idxPoint);
                end
            end
            alpha = 10.^(-alphaDB/20); % <-- original
%             alpha = 10.^(-0.05*alphaDB); % <-- original
%             alpha = exp(- 0.1151*alphaDB); % eq (1) from ISO 9613-1 - same res as above
            
%             alpha = 10.^(-alphaDB/10);
%             alpha = exp(-alphaDB/2);          
        end
        
        function R = TotalReflectionFactor(obj, eigenrays)
            %Calculates the total reflection factor of given eigenrays.
            %
            %This complex-valued factor is a catenation of the reflection
            %factors of successive reflections. If no reflection occured
            %the factor will be 1. The output has a size of nFreqs x nRays.
            %
            %Output:
            %R:         Complex valued factor containing information on
            %           phase and magnitude change of all reflections that
            %           occured for the given eigenray
            
            R = ones(numel(obj.frequencyVector), numel(eigenrays));
            for idxRay = 1:numel(eigenrays)
                eigenray = eigenrays(idxRay);
                R(:, idxRay) = obj.groundReflectionFactor.^eigenray.numReflections;                
            end
        end
    end
    
    %% Transfer Function (TF)
    methods( Access = public )
        function [combinedTF, individualTF] = TransferFunction(obj, eigenrays)
            %Generates an atmospheric transfer function based on given eigenrays.
            % Returns the joint TF of each eigenray as well as the combined
            % TF using the itaResult format.
            
            spreadingLoss = obj.SpreadingLoss(eigenrays);
            linearPhaseFilter = obj.PropagationDelayFilter(eigenrays);
            alpha = obj.AirAttenuation(eigenrays);
            R = obj.TotalReflectionFactor(eigenrays);
            
            individualTF = itaResult( alpha.*linearPhaseFilter.*R.*spreadingLoss, obj.frequencyVector, 'freq' );

            for idx = 1:numel(eigenrays)
                individualTF.channelNames{idx} = ['TF eigenray ' num2str(idx-1)];
            end
            
            combinedTF = itaResult( sum(individualTF.freqData, 2), obj.frequencyVector, 'freq' );
            combinedTF.channelNames = {'Atmospheric transfer function'};
        end
    end
    
    methods (Access = private)
        function linearPhaseFilter = PropagationDelayFilter(obj, eigenrays)
            delay = obj.PropagationDelay(eigenrays);
            linearPhaseFilter = exp(1j*2*pi.*delay.*obj.frequencyVector);
        end
    end
    
end