function [ absorption_db, alpha_iso ] = ita_atmospheric_absorption_level_dB( f, d, varargin )
	% ISO 9613-1 Acoustics - Attenuation of sound during propagation outdoors
        
	% Attenuation coefficient [dB/m], ~f, as assembly of Equation (5) parts
    if nargin == 2
        alpha_iso = ita_atmospheric_absorption_iso9613( f );
    else
        alpha_iso = ita_atmospheric_absorption_iso9613( f, varargin );
    end

	% Resulting atmospheric absorption [dB], ~alpha (~f) 
	% Equation (2)
    % Attenuation factor in decibel
	absorption_db = alpha_iso .* d;
    
end
