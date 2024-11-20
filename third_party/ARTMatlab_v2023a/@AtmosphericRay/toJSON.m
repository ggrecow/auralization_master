function jsonStr = toJSON(rays)
%AtmosphericRay.toJSON

% -------------------------------------------------------------------------
%                        ____  __________  _______
%                       //  / //__   ___/ //  _   |
%                      //  /    //  /    //  /_|  |
%                     //  /    //  /    //  ___   |
%                    //__/    //__/    //__/   |__|
%
% -------------------------------------------------------------------------
%                  ARTMatlab - ITAGeometricalAcoustics
%        (c) Copyright Institute of Technical Acoustics (ITA)
%  This file is part of the ARTMatlab application. Some rights reserved.
% You can find the license in the LICENSE.md file in the ARTMatlab folder.
%--------------------------------------------------------------------------

if numel(rays) == 0 || isempty(rays)
    warning('Cannot parse empty ray.')
elseif numel(rays) == 1
    jsonStruct = singleRayToJSONStruct(rays);
else
    jsonStruct = multipleRaysToJSONStruct(rays);
end

jsonStr = jsonencode(jsonStruct);

function jsonStruct = multipleRaysToJSONStruct(rays)
propagationPaths = cell(1, numel(rays));
for idxRay = 1:numel(rays)
    propagationPaths{idxRay} = singleRayToJSONStruct(rays(idxRay));
end

jsonStruct.class = 'propagation_path_list';
jsonStruct.identifier = '';
jsonStruct.propagation_paths = propagationPaths;

function jsonStruct = singleRayToJSONStruct(ray)

jsonStruct.class = 'propagation_path';
jsonStruct.identifier = '';

nPoints = ray.numPoints;
anchors = cell(1, nPoints);

anchors{1} = sourceJSONStruct(ray);
idxReflection = 1;
for idx = 2:nPoints-1
    isIdxReflection = any(ismember(ray.idxReflection, idx));
    if isIdxReflection
        anchors{idx} = specularReflectionJSONStruct(ray, idxReflection);
        idxReflection = idxReflection+1;
    else
        anchors{idx} = inhomogeneityJSONStruct(ray, idx);
    end
end
anchors{end} = receiverJSONStruct(ray);

jsonStruct.propagation_anchors = anchors;

function jsonStruct = initAnchorJSONStruct(ray, idxPoint)
jsonStruct.class = 'propagation_anchor';
r = ray.r.cart(idxPoint, :);
n = ray.n.cart(idxPoint, :);
if r(3) < 0
    r(3) = -r(3);
    n(3) = -n(3);
end
jsonStruct.interaction_point = [r'; 1];
jsonStruct.wavefront_normal = [n'; 1];
jsonStruct.time_stamp = ray.t(idxPoint);

function jsonStruct = inhomogeneityJSONStruct(ray, idxPoint)
jsonStruct = initAnchorJSONStruct(ray, idxPoint);
jsonStruct.anchor_type = 'inhomogeneity';

function jsonStruct = sourceJSONStruct(ray)
jsonStruct = initAnchorJSONStruct(ray, 1);
jsonStruct.anchor_type = 'inhomogeneity_source';

function jsonStruct = receiverJSONStruct(ray)
jsonStruct = initAnchorJSONStruct(ray, ray.numPoints);
jsonStruct.anchor_type = 'inhomogeneity_receiver';
if ~isempty(ray.mSpreadingLoss) && ray.mSpreadingLoss > 0
    jsonStruct.spreading_loss = ray.mSpreadingLoss;
end
if ~isempty(ray.isEigenray)
    jsonStruct.is_eigenray = ray.isEigenray;
end
if ~isempty(ray.receiverSphereHit)
    jsonStruct.receiver_hit = ray.receiverSphereHit;
end
if ~isempty(ray.rayZoomingIterations)
    jsonStruct.ray_zooming_iterations = ray.rayZoomingIterations;
end

function jsonStruct = specularReflectionJSONStruct(ray, idxReflection)
jsonStruct.class = 'propagation_anchor';
[r, n, t] = ray.reflectionParameters(idxReflection);
if r.z < 0
    r.z = -r.z;
    n.z = -n.z;
end
jsonStruct.interaction_point = [r.cart'; 1];
jsonStruct.wavefront_normal = [n.cart'; 1];
jsonStruct.time_stamp = t;

jsonStruct.face_normal = [0;0;1;1];
jsonStruct.anchor_type = 'inhomogeneity_specular_reflection';