function [resSpreadingLoss] = correctSpreadingLoss(atmoSpreadingLoss,dirPathLength,urbanPathLengths,isDiffracted,diffraction_model)
%CORRECTSPREADINGLOSS applies correction factor to atmospheric spreading 
% loss based on path length difference
% Syntax: [resSpreadingLoss] = correctSpreadingLoss(atmoSpreadingLoss,dirPathLength,urbanPathLengths,isDiffracted,diffraction_model)

if nargin < 5
    diffraction_model = 'utd';
end

%% calculate urban spreading loss
diffractedIDs = find(isDiffracted);

sumUrbanLength = sum(urbanPathLengths);




%% calculate correction factor based on path length difference
if isempty(diffractedIDs) || ~strcmp(diffraction_model,'utd')%if no diffraction during path
    
    correctionFactor = dirPathLength/sumUrbanLength;    
    urbanSpreadingLoss = 1/sumUrbanLength;
    
else%if diffracted path
    distToDiff = sum(urbanPathLengths(1:(diffractedIDs(1)-1)));
    distAfterDiff = sum(urbanPathLengths(diffractedIDs(1):length(urbanPathLengths)));
    

    urbanSpreadingLoss = 1/sqrt(distToDiff*distAfterDiff*(distToDiff+distAfterDiff));

    
    correctionFactor = dirPathLength/distToDiff * ...
         sqrt(distToDiff/(distAfterDiff*(distToDiff+distAfterDiff)));
    
    
%     if distToDiff > dirPathLength
%         correctionFactor = dirPathLength/distToDiff * ...
%             sqrt(distToDiff/(distAfterDiff*(distToDiff+distAfterDiff)));
%     else
%         correctionFactor = sqrt((dirPathLength-distToDiff)*dirPathLength/...
%             (distAfterDiff*(distToDiff+distAfterDiff)));
%     end
end

% caclulate correction factor between atmospheric and urban spreading loss
resSpreadingLoss = atmoSpreadingLoss/urbanSpreadingLoss * correctionFactor;

end

