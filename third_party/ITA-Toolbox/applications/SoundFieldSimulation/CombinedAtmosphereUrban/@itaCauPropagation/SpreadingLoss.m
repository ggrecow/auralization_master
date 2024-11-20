function [spreading_loss] = SpreadingLoss(obj, atmoSpreadingLoss, urbanPaths, diffraction_model)
%SPREADINGLOSS: Calculates spreading loss of a set of sound paths of a
%combined atmospheric urban simulation.
% All input variables are optional.

if nargin < 2; atmoSpreadingLoss = obj.calcAtmosphericParameters(); end
if nargin < 3; urbanPaths = obj.urbanPaths; end
if nargin < 4; diffraction_model = obj.diffractionModel; end

%% Derive parameters from paths
virtSourceDistance = norm( urbanPaths(1).propagation_anchors{end}.interaction_point -...
                           urbanPaths(1).propagation_anchors{1}.interaction_point );

%% Spreading loss
numPaths = numel(urbanPaths);
spreading_loss = zeros(1, numPaths);
for idPath = 1:numPaths
    
    nAnchors = numel(urbanPaths(idPath).propagation_anchors);
    isDiffracted = zeros(1,nAnchors-1);
    pathLengths = zeros(1,nAnchors-1);
    for idAnchor = 2:nAnchors
        
        pathLengths(idAnchor) = norm(urbanPaths(idPath).propagation_anchors{idAnchor}.interaction_point...
            - urbanPaths(idPath).propagation_anchors{idAnchor-1}.interaction_point);
        
        anchor_type = urbanPaths(idPath).propagation_anchors{idAnchor-1}.anchor_type;
        if strcmp(anchor_type,'outer_edge_diffraction') || strcmp(anchor_type,'inner_edge_diffraction')
            isDiffracted(idAnchor) = 1;
        end
        
    end
    
    spreading_loss(idPath) = adjust_spreading_loss(atmoSpreadingLoss,...
                                virtSourceDistance, pathLengths, isDiffracted, diffraction_model);
    
end

%% Adjust spreading loss
function [resSpreadingLoss] = adjust_spreading_loss(atmoSpreadingLoss,...
    dirPathLength, urbanPathLengths, isDiffracted, diffraction_model)
%Adjusts urban spreading loss based on atmospheric path

if nargin < 5
    diffraction_model = 'utd';
end

diffractedIDs = find(isDiffracted);
sumUrbanLength = sum(urbanPathLengths);

if isempty(diffractedIDs) || ~strcmp(diffraction_model,'utd')       
    correctionFactor = dirPathLength/sumUrbanLength;  
else   
    distToFirstDiff = sum(urbanPathLengths(1:(diffractedIDs(1)-1)));
    diffSpreadingLoss = 1/distToFirstDiff;
    lastInteraction = 1;
    for idDiff = 1:length(diffractedIDs)
                
        if (idDiff + 1) <= length(diffractedIDs)
            %consider segment to next aperture point
            nextAperture = diffractedIDs(idDiff+1);
        else
            %consider whole path after diffraction
            nextAperture = length(urbanPathLengths);
        end
        distToDiff = sum(urbanPathLengths(lastInteraction:(diffractedIDs(idDiff)-1)));
        distAfterDiff = sum(urbanPathLengths(diffractedIDs(idDiff):nextAperture));
        
        % make sure distances are ~= 0, assume geo precision
        distToDiff = max(distToDiff,1e-6);
        distAfterDiff = max(distAfterDiff,1e-6);
               
        diffSpreadingLoss = diffSpreadingLoss * sqrt(distToDiff/(distAfterDiff*(distToDiff+distAfterDiff)));
        lastInteraction = diffractedIDs(idDiff);
        
    end   
    correctionFactor = dirPathLength * diffSpreadingLoss;    
end

resSpreadingLoss = atmoSpreadingLoss * correctionFactor;
