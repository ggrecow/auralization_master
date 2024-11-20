classdef itaGeoPropagation < handle
    %ITAGEOPROPAGATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %---General properties---
        
        pps;                                        %Propagation paths struct
        c = 341.0;                                  %Speed of sound
        diffraction_model = 'utd';                  %Model for diffraction filter
        sim_prop = struct();                        %Struct with simulation settings like maximum reflection/diffraction order
        
        %---Transfer function properties---
        
        fs = 44100;                                 %Sampling frequency used for TF frequency vector
        
        %---Properties used to create auralization data---
        
        pps_old;                                    %Propagation paths of last time frame (only used for auralization parameters)
        source_id = 1;                              %Auralization source ID
        receiver_id = 1;                            %Auralization receiver ID
        freq_vector = ita_ANSI_center_frequencies'; %Frequency vector used for calculating auralization parameters
    end
    
   properties (Access = protected)
        n = 2^15 + 1;
        directivity_db = struct();
   end
    
    properties (Dependent)
        freq_vec;   %Linear spaced frequency vector used for TF calculation
        num_bins;   %Number of bins used for TF calculation
    end
    
    methods
        
        function obj = itaGeoPropagation( fs, num_bins )            
           if nargin >= 1
                obj.fs = fs;            
           end           
           if nargin >= 2
                obj.n = num_bins;            
           end
           
           obj.sim_prop.diffraction_enabled = true;
           obj.sim_prop.reflection_enabled = true;
           obj.sim_prop.directivity_enabled = true;
           obj.sim_prop.orders.reflection = -1;
           obj.sim_prop.orders.diffraction = -1;
           obj.sim_prop.orders.combined = -1;
           
        end
        
        function num_bins = get.num_bins( obj )
            num_bins = obj.n;
        end
        
        function f = get.freq_vec( obj )
            % Returns frequency base vector
            
            % taken from itaAudio (ITA-Toolbox)
            if rem( obj.n, 2 ) == 0
                f = linspace( 0, obj.fs / 2, obj.n )';
            else
                f = linspace( 0, obj.fs / 2 * ( 1 - 1 / ( 2 * obj.n - 1 ) ), obj.n )'; 
            end
            
        end
        
    end
end
