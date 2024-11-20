function [air_absorption] = AirAbsorption(obj, freqVector, atmoAirAbsorption, urbanPaths, constTemp, constRelHumidity, constStatPressure)
%AIRABSORPTION  Calculates corrected delay of a set of sound paths of a
%combined atmospheric urban simulation.
% All input variables are optional. Per default uses third-octave band
% center frequencies.


if nargin < 2; freqVector =  ita_ANSI_center_frequencies([20 20000], 3)'; end
if nargin < 3; atmoAirAbsorption = obj.calcAtmosphericParameters(freqVector); end
if nargin < 4; urbanPaths = obj.urbanPaths; end
if nargin < 5; constTemp = obj.constTemp; end
if nargin < 6; constRelHumidity = obj.constHumidity; end
if nargin < 7; constStatPressure = obj.constPressure; end

if isrow(atmoAirAbsorption)
    atmoAirAbsorption = atmoAirAbsorption.';
end
if isrow(freqVector)
    freqVector = freqVector';
end

%% Derive parameters from paths
virtSourceDistance = norm( urbanPaths(1).propagation_anchors{end}.interaction_point -...
                           urbanPaths(1).propagation_anchors{1}.interaction_point );
urbanPathLengths = ita_propagation_path_length( urbanPaths );
pathLengthDiff = urbanPathLengths - virtSourceDistance; %Difference compared to direct path

%% Calcualation
air_absorption = zeros(numel(freqVector), numel(urbanPaths));
for idPath = 1:numel(urbanPaths)
    air_absorption(:, idPath) = adjust_air_absorption(atmoAirAbsorption, pathLengthDiff(idPath), constTemp, constRelHumidity, constStatPressure, freqVector);
end

%% Adjust air absorption
function [outAbsorption] = adjust_air_absorption(inAbsorption, lengthDiff, T, humidity, p0, freqVector)
%Aplies correction factor to atmospheric absorption based on path length difference

absorp_factor = ita_atmospheric_absorption_iso9613(freqVector,T-273.15,humidity,p0/1000);
delta_dB = absorp_factor .* lengthDiff;

deltaAbsorption =  10.^( -delta_dB / 20.0 );

outAbsorption = inAbsorption .* deltaAbsorption;
