function receiverWFNormals = ReceiverWFNormals(urbanPaths)
%RECEIVERWFNORMALS Calculates the wavefront normals at the receiver for
%given urban paths (NX3 matrix).

%% Get WF normals
nPaths = numel(urbanPaths);
receiverWFNormals = zeros(nPaths, 3);
for idPath = 1:nPaths
    pathWFNormal = urbanPaths(idPath).propagation_anchors{end-1}.interaction_point(1:3)' -...
        urbanPaths(idPath).propagation_anchors{end}.interaction_point(1:3)';
    receiverWFNormals(idPath, :) = pathWFNormal / norm(pathWFNormal);
end
