function [transferFunction,transferFunction_separated] = TransferFunction(obj)
%TRANSFERFUNCTION Calculates combined transfer function based on an
%atmospheric freefield path and a set of urban paths calculated with the
%virtual source method.
%   Note: The members atmosphericPath and urbanPaths have to be set first.

%% Check for sound paths and atmosphere settings
assert(~isempty(obj.atmosphericPath) && isa(obj.atmosphericPath, 'AtmosphericRay'), '[itaCauPropagation]: Atmospheric path is not set properly')
assert(~isempty(obj.urbanPaths) && isstruct(obj.urbanPaths), '[itaCauPropagation]: Urban paths are not set')
assert( ~isempty(obj.atmosphere) && isa(obj.atmosphere, 'StratifiedAtmosphere'), '[itaCauPropagation]: Stratified atmosphere is not set properly.' )

%% Read lengths and delay of urban Paths
numPoints = numel(obj.urbanPaths(1).propagation_anchors);

virtSource = obj.urbanPaths(1).propagation_anchors{1}.interaction_point;
receiverPos = obj.urbanPaths(1).propagation_anchors{numPoints}.interaction_point;
virtSourceDistance = norm(virtSource - receiverPos);

[pathLengths,timeDelay,isDiffracted] = extractUrbanPathLengths(obj.urbanPaths);

%% Make sure nBins is sufficient for path lengths
fftDegree = ceil(log2(max(sum(timeDelay,2))*obj.samplingRate));
obj.nBins = max(obj.nBins,2^fftDegree);
freqVector = obj.getFreqVector;

%% Calculate urban TF
c = obj.constSpeed;
numPaths = numel(obj.urbanPaths);
[~,urbanTF_separated] = obj.calcUrbanTF;

%% Atmospheric TF parameters
[spreading_loss, propagation_delay, air_attenuation] = obj.calcAtmosphericParameters(freqVector);

%% Adapt atmospheric parameters based on difference in path length and delay
adaptedSpreadingLoss = ones(1,numPaths);
deltaPhase = ones(numPaths,length(freqVector));
adaptedAbsorption = ones(numPaths,length(freqVector));

for idPath=1:numPaths

    if strcmp(obj.virtualSourceModus,'delay')
        [adaptedSpreadingLoss(idPath)] = correctSpreadingLoss(spreading_loss,virtSourceDistance,...
                  pathLengths(idPath,:),isDiffracted(idPath,:), obj.diffractionModel);
    else 
        deltaPhase(idPath,:) = correctPhaseDelay(propagation_delay,virtSourceDistance,freqVector,c);
    end

    adaptedAbsorption(idPath,:) = obj.AirAbsorption( freqVector, air_attenuation,...
        obj.urbanPaths(idPath), obj.constTemp, obj.constHumidity,obj.constPressure );

end

% Combine attenuation from air absorption and spreading loss and phase delay due to propagation
atmosTF_separated = adaptedSpreadingLoss' .* adaptedAbsorption .* deltaPhase;

%% Apply atmospheric parameters to urban transfer function
[transferFunction,transferFunction_separated] = applyAtmoAttenuation(urbanTF_separated,atmosTF_separated,freqVector);

transferFunction.channelNames = {'Combined TF'};

for idPath = 1:transferFunction_separated.nChannels
   transferFunction_separated.channelNames{idPath} = ['TF for path ',num2str(idPath)]; 
end
