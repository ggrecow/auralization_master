classdef itaCauPropagation < handle
%ITACAUPropagation Calculates auralization parameters for combined
%propagation model
    
    properties(Access = public)
        % ---Sound paths---
        
        atmosphericPath = [];   %Atmospheric path used for calculations
        urbanPaths = [];        %Urban sound paths used for calculations
        
        % ---Inhomogeneous medium---
        
        atmosphere = []; %StratifiedAtmosphere object from ARTMatlab
        
        % ---Properties for constant medium---
        
        constSpeed = 340;       %[m/s]
        constTemp = 20;         %[°C]
        constHumidity = 50;     %Relative humidity [%]
        constPressure = 101325; %Static pressure [Pa]
        
        % ---Settings---
        
        virtualSourceModus = 'delay';
        maxReflectionOrder = 1;
        
        diffractionModel = 'utd';
        nBins = 2^20;
        
    end
    
    properties(Access = public, Hidden = true)
        
        samplingRate = 44100;
        
    end
    
    %% Initialization
    methods
        function obj = itaCauPropagation( cauSimulationEngine )
            %Creates an object of the CAU propagation model. Optionally,
            %you can hand a CAU simulation engine to transfer the simulation
            %settings (see TransferSettings()).
            if nargin
                assert( isa(cauSimulationEngine, 'itaCauSimulation'), 'Input must be an itaCauSimulation object.' );
                obj.TransferSettings( cauSimulationEngine );
            end
        end
   end 
    %% public functions
    methods(Access = public)
        
        function TransferSettings(obj, cauSimulationEngine)
            %Transfers simulation settings from an itaCauSimulation object
            %including StratifiedAtmosphere
            
            assert( isa(cauSimulationEngine, 'itaCauSimulation'), 'Input must be an itaCauSimulation object.' );
            
            %Simulation settings
            obj.virtualSourceModus = cauSimulationEngine.virtualSourceModus;
            obj.maxReflectionOrder = cauSimulationEngine.urbanPropagation.MaxReflectionOrder;
            
            %Stratified Atmosphere
            obj.atmosphere = cauSimulationEngine.atmosphere;
        end
        function SetSoundPaths(obj, atmoPath, urbanPaths, c, hRel, T, p0)
            %Sets the sound paths from a CAU simulation which can then be
            %used to calculate auralization parameters or a transfer
            %function.
            %Optionally, the parameters for the homogeneous medium can be
            %set. Per default, the values are taken from the stratified
            %atmosphere at receiver altitude.
            assert( isa(atmoPath, 'AtmosphericRay'), 'First input must be of type AtmosphericRay.' );
            assert( isstruct(urbanPaths) && isfield(urbanPaths, 'propagation_anchors'), 'Second input must be a set of urban sound paths.' );
            
            %Sound paths
            obj.atmosphericPath = atmoPath;
            obj.urbanPaths = urbanPaths;
            
            %---Homogeneous atmosphere---
            %Default: Use parameters at receiver altitude
            receiverPosition = urbanPaths(1).propagation_anchors{end}.interaction_point;
            receiverAltitude = receiverPosition(3);
            
            if nargin < 4
                obj.constSpeed = obj.atmosphere.c(receiverAltitude);
            else
                obj.constSpeed = c;
            end
            if nargin < 5
                obj.constHumidity = obj.atmosphere.humidity(receiverAltitude);
            else
                obj.constHumidity = hRel;
            end
            if nargin < 6
                obj.constTemp = obj.atmosphere.T(receiverAltitude);
            else
                obj.constTemp = T;
            end
            if nargin < 7
                obj.constPressure = obj.atmosphere.p0(receiverAltitude);
            else
                obj.constPressure = p0;
            end
            
        end
        
        
        [combinedTF_sum,combinedTF_separated] = TransferFunction(this);
        
        function [freqVector] = getFreqVector(obj)
            
            if rem( obj.nBins, 2 ) == 0
                freqVector = linspace( 0, obj.samplingRate / 2, obj.nBins )';
            else
                freqVector = linspace( 0, obj.samplingRate / 2 * ( 1 - 1 / ( 2 * obj.nBins - 1 ) ), obj.nBins )';
            end
            
        end
        
        [airAbsorption, objectInteraction,spreadingLoss,diffractionGain,delay,sourceWFNormal,receiverWFNormal] = ...
            AuralizationParameters(obj, freqVector, atmoPropDelay)
        [auralization_frame_files] = ...
            AuralizationScenarioParameters(obj, sim_result_files, outputFolder, freqVector, force)
        
        [spreading_loss] = SpreadingLoss(obj, atmoSpreadingLoss, urbanPaths, diffraction_model)
        [delay] = PropagationDelay(obj, atmoDelay, urbanPaths, constSpeedOfSound)
        [air_absorption] = AirAbsorption(obj, freqVector, atmoAbsorption, urbanPaths, constTemp, constRelHumidity, constStatPressure)
    end
    
    methods(Access = private)
        function [spreading_loss, propagation_delay, air_attenuation] = calcAtmosphericParameters(obj, frequencyVector)
            
            assert( ~isempty(obj.atmosphericPath) && isa(obj.atmosphericPath, 'AtmosphericRay'), '[itaCauPropagation]: Atmospheric path is not set properly' )
            assert( logical(exist('AtmosphericPropagation', 'class')), '[itaCauPropagation]: Cannot find AtmosphericPropagation class from ARTMatlab. Did you make it available in Matlab path?' )
            assert( ~isempty(obj.atmosphere) && isa(obj.atmosphere, 'StratifiedAtmosphere'), '[itaCauPropagation]: Stratified atmosphere is not set properly.' )
            
            spreading_loss = obj.atmosphericPath.spreadingLoss;
            propagation_delay = obj.atmosphericPath.t(end);
            
            if nargout > 2
                atmosPropModel = AtmosphericPropagation;
                atmosPropModel.atmosphere = obj.atmosphere;
                atmosPropModel.frequencyVector = frequencyVector;
                air_attenuation = atmosPropModel.AirAttenuation(obj.atmosphericPath);
            end
        end
        
        function [urbanTF_sum,urbanTF_separated] = calcUrbanTF(obj)
        
            % initialize itaGeoPropagation 
            numBins = obj.nBins;
            itaInfiniteWedge.set_get_geo_eps(1e-6);

            gpsim = itaGeoPropagation(obj.samplingRate,numBins);
            gpsim.c = obj.constSpeed;

            gpsim.diffraction_model = obj.diffractionModel;

            urbanTF_sum = itaAudio;
            urbanTF_separated = itaAudio;

            % calculate separate and combined transfer function
            for idx=1:numel(obj.urbanPaths)

                gpsim.pps = obj.urbanPaths(idx);

                pathTF = itaAudio();
                try
                    pathTF.freqData = gpsim.run;
                catch err
                    warning('Could not determine transfer function for path %i:"%s"', idx, err.message)
                    pathTF.freqData = 1e-12*ones(numBins,1);
                end
                if idx == 1
                    urbanTF_separated = pathTF;
                    urbanTF_sum = pathTF;
                else
                    urbanTF_separated = ita_merge( urbanTF_separated, pathTF );
                    urbanTF_sum = ita_add( urbanTF_sum, pathTF );
                end

            end
        end
    end
    
    
end

