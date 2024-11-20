function ray = parseJSON(jsonStr)
%AtmosphericRay.parseJSON

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

%% Decode
jsonStruct = jsondecode(jsonStr);

%% Parse individual profiles
switch(lower(jsonStruct.class))
    case 'propagation_path'
        ray = parseToSingleRay(jsonStruct);
    case 'propagation_path_list'
        ray = parseToMultipleRays(jsonStruct);
    otherwise
        error('Class "%s" is not supported.', jsontStruct.class)
end

function rays = parseToMultipleRays(jsonStruct)

rays = [];
for idxRay = 1:numel(jsonStruct.propagation_paths)
    pathStruct = jsonStruct.propagation_paths(idxRay);
    rays = [rays parseToSingleRay(pathStruct)];
end

function ray = parseToSingleRay(jsonStruct)

ray = AtmosphericRay();
for idxAnchor = 1:numel(jsonStruct.propagation_anchors)
    if iscell(jsonStruct.propagation_anchors)
        anchorStruct = jsonStruct.propagation_anchors{idxAnchor};
    else
        anchorStruct = jsonStruct.propagation_anchors(idxAnchor);
    end
    assert(strcmpi(anchorStruct.class, 'propagation_anchor'), 'Expected an element of class propagation_anchor instead of "%s"', anchorStruct.class)
    assert(contains(lower(anchorStruct.anchor_type), 'inhomogeneity'), 'Expected an anchor with any "inhomogeneity" type. Instead found type "%s"', anchorStruct.anchor_type)
    
    r = anchorStruct.interaction_point(1:3)';
    n = anchorStruct.wavefront_normal(1:3)';
    t = anchorStruct.time_stamp;    
    switch lower(anchorStruct.anchor_type)
        case 'inhomogeneity_source'
            if abs(r(3)) < eps && n(3) < 0 %Source at z=0 and shooting below surface
                ray.addReflection(r, n, 0);
            else
                ray.init(r,n);
            end
        case 'inhomogeneity_receiver'
            ray.addData(r, n, t);
            if isfield(anchorStruct, 'spreading_loss') && anchorStruct.spreading_loss > 0
                ray.mSpreadingLoss = anchorStruct.spreading_loss;
            elseif isfield(anchorStruct, 'spreading_loss_db')
                ray.mSpreadingLoss = 10^(anchorStruct.spreading_loss_db/20);
            end
            if isfield(anchorStruct, 'is_eigenray')
                ray.isEigenray = anchorStruct.is_eigenray;
            end
            if isfield(anchorStruct, 'receiver_hit')
                ray.receiverSphereHit = anchorStruct.receiver_hit;
            end
            if isfield(anchorStruct, 'ray_zooming_iterations')
                ray.rayZoomingIterations = anchorStruct.ray_zooming_iterations;
            end
        case 'inhomogeneity'
            ray.addData(r, n, t);
        case 'inhomogeneity_specular_reflection'
            ray.addReflection(r, n, t);
        otherwise
            error('Anchor type "%s" is not supported for parsing into AtmosphericRay', anchorStruct.anchor_type)
    end
end