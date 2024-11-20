function sourceWFNormals = SourceWFNormals(urbanPaths,atmoPath)
%SOURCEWFNORMALS Calculates the wavefront normals at the source for given
%urban paths. The wavefront normals are rotated according to
%the difference between the source WF normal of the atmospheric path and
%the direct urban path.

%% Rotation between source WF normal of real and virtual direct path
virtualSource = urbanPaths(1).propagation_anchors{1}.interaction_point(1:3)';
receiver = urbanPaths(1).propagation_anchors{end}.interaction_point(1:3)';
virtualWFNormal = receiver - virtualSource;
virtualWFNormal = virtualWFNormal / norm(virtualWFNormal);

atmoWFNormal = atmoPath.n0.cart;
atmoWFNormal = atmoWFNormal / norm(atmoWFNormal);

rotMat = RotationMatrix(virtualWFNormal, atmoWFNormal);

%% Get WF normals and rotate them
nPaths = numel(urbanPaths);
sourceWFNormals = zeros(nPaths, 3);
for idPath = 1:nPaths
    nAnchors = numel(urbanPaths(idPath).propagation_anchors);
    if nAnchors == 2 %Direct path => take original WF normal
        sourceWFNormals(idPath, :) = atmoWFNormal;
    else %Other paths
        pathWFNormal = urbanPaths(idPath).propagation_anchors{2}.interaction_point(1:3)' -...
            urbanPaths(idPath).propagation_anchors{1}.interaction_point(1:3)';
        pathWFNormal = pathWFNormal / norm(pathWFNormal);
        
        sourceWFNormals(idPath, :) = pathWFNormal * rotMat;
    end
end

function R = RotationMatrix(virtualWFNormal, atmoWFNormal)

C = cross(virtualWFNormal, atmoWFNormal) ;
D = dot(virtualWFNormal, atmoWFNormal) ;
norm_factor = norm(virtualWFNormal) ; % used for scaling

if ~all(C==0) % check for colinearity
    Z = [0 -C(3) C(2); C(3) 0 -C(1); -C(2) C(1) 0] ;
    R = (eye(3) + Z + Z^2 * (1-D)/(norm(C)^2)) / norm_factor^2 ; % rotation matrix
else
    R = sign(D) * (norm(atmoWFNormal) / norm_factor) ; % orientation and scaling
end

R = R'; % Since we are using row vectors
