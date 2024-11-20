function [weights] = ita_pan_vbap(LS_pos, VS_pos, varargin)
%ITA_PAN_VBAP Use the VBAP method from
% References:
%   [1] Pulkki; Virtual sound source positioning using vactor base amplitude panning; 1997
%
% Arguments:
%   LS_pos: itaCoordinates, loudspeaker setup
%   VS_pos: itaCoordinates, virtual source positions
% Options: (default)
%   {'method'      : ('static') , use convex hull for triangulation
%                    'dynamic'  , calculate loudspeaker triangle based on distance
%    'normalize'   : (true)     , normalize weigths
%                    false      , NOTE: Distance Loss will not be applied if normalize is false
%    'distanceLoss': (true)     , include 1/r - law
%                    false      , panning only}
%
% Returns:
%   weights: VBAP matrix

%% Parse input
params.method = 'static';
params.normalize = true;
params.distanceLoss = true;
[params, ~] = ita_parse_arguments(params, varargin);

%% Initialize
weights = zeros(VS_pos.nPoints, LS_pos.nPoints);
LS_aux = LS_pos;
LS_aux.r = 1;
VS_aux = VS_pos;
VS_aux.r = 1;

%% Calculate weights
% Static Approach
if strcmp(params.method, 'static')
    % calculate triangles
    H = convhull(LS_pos.x, LS_pos.y, LS_pos.z);
    % unhandled virtual sources
    unhandled_VS = true(VS_aux.nPoints, 1);
    % calculate weights for each source
    for hdx = 1:size(H, 1)
        % stop loop if all VS are handled
        if sum(unhandled_VS) == 0
            break
        end
        
        % calculate gains
        p = VS_aux.cart;
        L = LS_aux.n(H(hdx, :)).cart;
        g = p * pinv(L);
        
        % find correct triangle
        idx = all(g > -1e-5, 2) & unhandled_VS;
        g = g(idx, :);
        g(abs(g) <= 1e-5) = 0;
        % normalize
        if params.normalize && any(idx)
            activity = 1;  % ./ sqrt(sum(g > 0, 2));
            g = activity .* abs(g) ./ vecnorm(g, 2, 2);
            % account for difference in distance loss
            if params.distanceLoss
                g = (LS_pos.n(H(hdx, :)).r' ./ VS_pos.r(idx)) .* g;
            end
        end
        
        % write weights to matrix
        weights(idx, H(hdx, :)) = g;
        unhandled_VS(idx) = false;
    end

    
% Dynamic Approach
elseif strcmp(params.method, 'dynamic')
    warning('ita_pan_vbap: only use the dynamic approach if you know what it is doing. May result in non-optimal VBAP triangles.')
    for idx = 1:VS_pos.nPoints
        % get loudspeakers with shortest distance to line between origin and source
        t = LS_aux.dot(repmat(VS_aux.n(idx), LS_aux.nPoints));
        valid = find(t > -1e-15);
        x = itaCoordinates(LS_aux.cart - t.*VS_aux.cart(idx, :));
        % choose closest 3 points to line in the direction of the line (t<0)
        [~, jdx] = sort(x.r(valid), 'ascend');
        % select a non coplanar triplet as initial triplet
        triplets = nchoosek(jdx, 3);
        for tdx = 1:size(triplets, 1)
            jdx_active = valid(triplets(tdx, :));
            LS_active = LS_aux.n(jdx_active);
            if abs(LS_active.n(1).dot(LS_active.n(2).cross(LS_active.n(3)))) > 0.1
                break
            end
            % throw error when all possibilities are invalid
            if tdx == size(triplets, 1)
                error('All initial triplets are coplanar!')
            end
        end
        
        % test first selection
        p = VS_aux.cart(idx, :);
        g = p * pinv(LS_active.cart);
        % check for negative gains
        while_counter = 0;
        while any(g < -1e-15)
            if while_counter > LS_aux.nPoints
                error('Infinite loop detected')
            end
            % adapt points whose weights are negative
            for fdx = find(g < -1e-15)
                % get opposite direction distance vector from line
                x_new = -x.n(jdx_active(fdx));
                x_new.r = 1;
                % calculate base point of that vector
                p_ref = itaCoordinates(t(jdx_active(fdx))*p);
                % shift other loudspeaker points by base
                LS_shift = LS_aux - p_ref;
                % prefer close points with bad scalar product to far
                % sources with good scalar product
                LS_shift.r = 1./LS_shift.r;
                % use point that best matches direction as new point
                [~, mdx] = sort(LS_shift.dot(repmat(x_new, LS_shift.nPoints)), 'descend');
                % select a fitting mdx
                for m = 1:numel(mdx)
                    if ~any(jdx_active == mdx(m))
                        jdx_active(fdx) = mdx(m);
                        LS_active.cart(fdx, :) = LS_aux.cart(mdx(m), :);
                        if abs(LS_active.n(1).dot(LS_active.n(2).cross(LS_active.n(3)))) > 0.1
                            break
                        end
                    end
                    % Throw error when no fitting mdx is found
                    if m == numel(mdx)
                        error('All alternative LSs are coplanar!')
                    end
                end
                
                % check if only moving one point solves problem
                g = p * pinv(LS_active.cart);
                if all(g > -1e-15)
                    break
                end
            end
            % no need to update g, as it is updated in for loop
            while_counter = while_counter + 1;
        end
        g(abs(g) < 1e-5) = 0;
        
%         % get loudspeaker groups with equal distance
%         d_change = find(~(diff(d) < 0.03));
%         d_groups = diff([0, d_change']);
%         n_groups = find(cumsum(d_groups) >= 3, 1);
%         
%         % construct all possible loudspeaker triangles
%         
%         if n_groups == 1 % 3 or more with equal distance
%             tris = nchoosek(map(1:d_groups(1)), 3);
%         elseif n_groups == 2 % first or first two with equal distance
%             % construct triangles where group 1 is a fixed part of all
%             % triangles
%             selection = d_groups(1)+1 : d_groups(1)+d_groups(2);
%             k = abs(3 - d_groups(1));
%             tris = [repmat(map(1:d_groups(1))', nchoosek(d_groups(2), k), 1), ...
%                     nchoosek(map(selection), k)];
%         else % n_groups == 3, there should be no other posibilities
%             selection = 3 : sum(d_groups(1:3));
%             tris = [repmat(map(1:2).', numel(selection), 1), map(selection)];
%         end
%         
%         % calculate gains for all triangles
%         g = zeros(size(tris));
%         for tdx = 1:size(tris, 1)
%             jdx_active = tris(tdx, :);
%             LS_active = LS_aux.n(jdx_active);
%             if abs(LS_active.n(1).dot(LS_active.n(2).cross(LS_active.n(3)))) < 1e-15
%                 g(tdx, :) = nan;
%                 continue
%             end
%             
%             % calculate gains
%             p = VS_aux.cart(idx, :);
%             L = LS_active.cart;
%             g(tdx, :) = p * pinv(L);
%         end
%         
%         % select triangle based on stability first
%         % [~, tdx] = min(max(g, [], 2)) may be sharpness first
%         [~, tdx] = max(min(g, [], 2));
%         jdx_active = tris(tdx, :);
%         g = g(tdx, :);
%         g(abs(g) <= 1e-5) = 0;
%         LS_active = LS_pos.n(jdx_active);
        
        % normalize
        if params.normalize
            % distribute energy to all
            activity = 1/sqrt(sum(g > 0));
            g = activity * abs(g)/ norm(g);
            % account for difference in distance loss
            if params.distanceLoss
                g = LS_active.r' / VS_pos.r(idx) .* g;
            end
        end
        % write weights to matrix
        weights(idx, jdx_active) = g;
    end
else
    error("I do not know that method! Try 'static' or 'dynamic'.");
end

% transpose for easier usage
weights = weights.';

end

