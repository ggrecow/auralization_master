function [pathLengths,timeDelay,isDiffracted,virtSourceDistance] = extractUrbanPathLengths(urbanPaths,c)
%EXTRACTURBANPATHLENGTHS extracts lengths and delays of individual paths
%segments from urban Paths

if nargin < 2
    c = 340;
end

numPaths = numel(urbanPaths);
maxAnchors = 9;

pathLengths = zeros(numPaths,maxAnchors);
isDiffracted = zeros(numPaths,maxAnchors); %indicates if path segment is diffracted
virtSourceDistVector = zeros(numPaths,1);

for idPath = 1:numPaths
    
   for idAnchor = 2:numel(urbanPaths(idPath).propagation_anchors)
       
           pathLengths(idPath,idAnchor) = norm(urbanPaths(idPath).propagation_anchors{idAnchor}.interaction_point...
               - urbanPaths(idPath).propagation_anchors{idAnchor-1}.interaction_point);

           anchor_type = urbanPaths(idPath).propagation_anchors{idAnchor-1}.anchor_type;
           if strcmp(anchor_type,'outer_edge_diffraction') || strcmp(anchor_type,'inner_edge_diffraction')
                isDiffracted(idPath,idAnchor) = 1;
           end
       
   end
   virtSourceDistVector(idPath) = norm(urbanPaths(idPath).propagation_anchors{end}.interaction_point-...
       urbanPaths(idPath).propagation_anchors{1}.interaction_point);
   
end

virtSourceDistance = mean(virtSourceDistVector);
timeDelay = pathLengths/c;

end

