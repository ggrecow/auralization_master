function [airAbsorption, objectInteraction, spreadingLoss, diffractionGain, delay, sourceWFNormals, receiverWFNormals] = AuralizationParameters(obj, freqVector, atmoPropDelay)
%AURALIZTIONPARAMETERS Calculates parameters used for auralization
% For all results, the second dimension refers to the index of the
% respective sound paths. The first dimension, refers to frequency for
% frequency-dependent parameters.

if nargin < 2; freqVector = ita_ANSI_center_frequencies([20 20000], 3); end

%Ensure column vector for frequencies
if isrow(freqVector); freqVector = freqVector'; end

numPaths = numel(obj.urbanPaths);
nBins = numel(freqVector);

%% Atmospheric parameters
[atmoSpreadingLoss, atmoDelay, atmoAirAbsorption] = obj.calcAtmosphericParameters(freqVector);
if nargin >= 3
    atmoDelay = atmoPropDelay;
end

%% Final paramters: Spreading loss, delay, air absorption
spreadingLoss = obj.SpreadingLoss( atmoSpreadingLoss );
delay = obj.PropagationDelay( atmoDelay );
airAbsorption = obj.AirAbsorption( freqVector, atmoAirAbsorption );
    
%% Diffraction filter
diffractionFilter = zeros(nBins, numPaths);
diffractionGain = zeros(1, numPaths);
for idPath=1:numPaths
    [diffractionFilter(:,idPath), diffractionGain(idPath)] = calcDiffractionFilter(obj.urbanPaths(idPath), freqVector, obj.constSpeed, obj.maxReflectionOrder);
end

%% Reflection
%TODO...

%% Accumulated third-octave magnitudes caused by object interaction
objectInteraction = abs(diffractionFilter);
if any( objectInteraction == 0 )
    warning('Frequency magnitudes not determined correctly')
end

%% Correct gain factor based on diffraction
% For phase inversion due to diffraction (+1, -1) or if diffraction path is not compensated (0)
diffractionGain = sign(diffractionGain);

%% WF normals
if nargout >= 5
    sourceWFNormals = SourceWFNormals(obj.urbanPaths, obj.atmosphericPath);
end

if nargout >= 6
    receiverWFNormals = ReceiverWFNormals(obj.urbanPaths);
end
