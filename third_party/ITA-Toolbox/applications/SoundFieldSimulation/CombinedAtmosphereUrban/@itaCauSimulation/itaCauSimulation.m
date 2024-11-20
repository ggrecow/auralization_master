classdef itaCauSimulation
    %ITA_CAU_AURALIZATION interface class for the combined atmospheric and
    %urban auralization
    
    properties(Hidden = true)
        virtualSourceModus = 'delay';
    end
    properties
        windDirection = [-1 0 0];
        atmosphere
        rayTracer
        urbanPropagation
        cityModel
        filter_delays = true;
    end
    
    %% initialization
    methods
        function obj = itaCauSimulation()
            %ITA_CAU_AURALIZATION Construct an instance of this class     
            obj.rayTracer = initializeRayTracer;
            obj.atmosphere = initializeAtmosphere([-1 0 0]);          
            obj.urbanPropagation = initializeUrbanPropagation;
        end
        
    end
    
    %% Hidden functions
    methods(Access = public, Hidden = true)
        
        function [atmoPaths, receiverMissed] = RunRayTracing(obj,sourcePosition,receiverPosition)
            %RUNRAYTRACING Find eigenray for given source and receiver position using the Atmospheric Ray Tracing framework
            receiverMissed = 1;
            atmoPaths = [];
            try
                atmoPaths = obj.rayTracer.FindEigenrays(obj.atmosphere,sourcePosition,receiverPosition);
            catch
                tmpAngle = obj.rayTracer.maxAngleForGeomSpreading;
                obj.rayTracer.maxAngleForGeomSpreading = 0.1;
                try 
                    atmoPaths = obj.rayTracer.FindEigenrays(obj.atmosphere,sourcePosition,receiverPosition);
                catch
                end
                obj.rayTracer.maxAngleForGeomSpreading = tmpAngle;
            end
            
            if isempty(atmoPaths)
                warning('WARNING: Error during calculation of eigenrays. No eigenray found.')
                return;
            else
                if isprop(atmoPaths(1), 'receiverSphereHit')
                    receiverMissed = ~atmoPaths(1).receiverSphereHit;
                else
                    minDistancePoint = atmoPaths(1).r.cart(end,:);
                    distToReceiver = norm(minDistancePoint-receiverPosition);
                    receiverMissed = distToReceiver > obj.rayTracer.maxReceiverRadius;
                end
                receiverMissed = double( receiverMissed );
            end
         
            if receiverMissed
                warning('WARNING: Receiver not reached by eigenray.')
            end
        end

        function [urbanPaths] = RunUrbanPropagation(obj,sourcePosition,receiverPosition)
            %RUNURBANPROPAGATION Calculation of urban sound propagation
            %paths for given source and receiver position using pigeon
            urbanPaths = obj.urbanPropagation.run(sourcePosition,receiverPosition);
        end
    end
    
    %% Public functions
    methods(Access = public)
        function [atmoPath, virtualSource, urbanPaths] = Run(obj,sourcePos,receiverPos)
        %RUN Runs combined simulation and returns the atmospheric free
        %field path, the virtual source positions and the resulting urban
        %sound paths

            % Run ray tracing
            [atmoPath, receiverMissed] = obj.RunRayTracing(sourcePos,receiverPos);

            if receiverMissed
                error('ART: Eigenrays could not be determined.')
            end

            virtualSource = obj.VirtualSource(atmoPath,receiverPos);
            
            % Run urban propagation
            if nargout > 2
                urbanPaths = obj.urbanPropagation.run(virtualSource,receiverPos);
            end
        end
        
        sim_result_files = DynamicScenario(obj, sourcePositions, receiverPosistions, updateRate, outputFolder, force)
    end
    
     %% Private functions
    methods(Access = public)
       
        [virtualSource] = VirtualSource(obj, atmoPath, receiverPosition, delay)
        
    end
    
    
    %% set functions
    methods
        
        function obj = set.windDirection(obj,input)
            if ~isnumeric(input) || ~isequal(size(input),[1 3])
                error('Wind direction must be defined as 1x3 numeric vector.')
            end
            obj.windDirection = input;
            obj.atmosphere = initializeAtmosphere(input);          
        end
        
        function obj = set.cityModel(obj,input)
            obj.cityModel = input;
            obj.urbanPropagation.geometry_file_path = input;         
        end
        
    end
    
end

