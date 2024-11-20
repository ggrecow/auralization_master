function [contours] = ita_setup_checkForHoles(LS_pos, varargin)
%ITA_SETUP_CHECKFORHLES Check the convex hull of the loudspeaker setup for holes
%   The convex hull is searched for holes, which are defined by a max 
%   aperture angles between two loudspeakers (from setup origin).
%
% Arguments:
%   LS_pos: itaCoordinates, loudspeaker setup
% Options: (default)
%   {H          : ([]), hull of loudspeaker setup given in index triplets
%    maxAperture: (90), triplets with an aperture >= this value are considered holes}
%
% Returns:
%   LS_imag: itaCoordinates, position of imaginary loudspeakers

%% Parse input
params.H = [];
params.maxAperture = 90;
params = ita_parse_arguments(params, varargin);

% check convex hull H
if isempty(params.H)
    H = convhull(LS_pos.x, LS_pos.y, LS_pos.z);
elseif size(params.H, 2) ~= 3
    warning('Given hull is not made of triplets, calculating own hull using convhull!')
    H = convhull(LS_pos.x, LS_pos.y, LS_pos.z);
else
    H = params.H;
end

%% Search for holes in convex hull
% store points, that lie on the contours of 'holes'
contours = {};
% check each trinagle
for idx = 1:size(H, 1)
    % calculate angle at the listener from distance between speakers
    tri = H(idx, :);
    d = [LS_pos.n(tri(1)) - LS_pos.n(tri(2)) ...
         LS_pos.n(tri(2)) - LS_pos.n(tri(3))...
         LS_pos.n(tri(3)) - LS_pos.n(tri(1))];
    d = d.merge;
    alpha = 2*asin(d.r/(2*LS_pos.r(1)));
    % check if any angle of a given loudspeaker triple is >= maxAperture
    if any(alpha >= params.maxAperture * pi/180)
        % start new contour, if empty
        if isempty(contours)
            contours = {tri};
        else % add triangle to an existing contour or star a new one
            for jdx = 1:numel(contours)
                % check if two of the triangle points are already on the contour
                % contours intersecting at one point will be viewed as two separate contours
                if numel(intersect(contours{jdx}, tri)) == 2
                    % Add missing points to contour
                    contours{jdx} = unique([contours{jdx}, tri]);
                    tri = [];
                    break;
                end
            end
            if ~isempty(tri) % no matching contour, start a new one
                contours{end+1} = tri;
            end
        end
    end
end

% merge contours, as a contour may have evolved from different starting points
if ~isempty(contours)
    contours = merge_contours(contours);
end

end

%% Merge contours helper function
function [merged_contours] = merge_contours(contours)
%MERGE_CONTOURS Takes contours as a cell of index arrays and merges these
%arrays if two or more indices are equal
%   Stepwise merging of the contours by comparing all posible combinations
%   of two contours. If the contours are mergable, they are merged and
%   written to a new cell element of merged contours. After all
%   combinations are done, the contours used for merging are removed from
%   the input cell. The unused contours are not mergable and will be added
%   to the output later. Now the merging is repeated on the new cell, until
%   no merges are possible or only one contour is left. At the end of each
%   recursion the unmerged contours get added back to the output.
%
%   Example: contours = {[1,2,3], [2,3,4,5,6], [5,6,7], [4, 8, 9]}
%       1. step: merged   = {[1,2,3,4,5], [2,3,4,5,6,7]}
%                unmerged = {[4,8,9]}
%       2. step: merged   = {[1,2,3,4,5,6,7]}
%                unmerged = {}
%       3. step: output   = {[1,2,3,4,5,6,7], [4,8,9]}
%
% Arguments:
%   contours: cell, contains arrays of point indices as its elements
%
% Returns:
%   merged_contours: cell, containes the merged contours

%% Recursively merge contours
% Only one contour, nothing to do
if numel(contours) == 1
    merged_contours = contours;
else
    % Initialize
    contour_pairs = nchoosek(contours, 2);
    idx = nchoosek(1:numel(contours), 2);
    merged_contours = {};
    merged_idx = [];
    n_merges = 0;
    % Merge a contour pair if two or more points are identical
    for jdx = 1:size(contour_pairs, 1)
        if numel(intersect(contour_pairs{jdx, 1}, contour_pairs{jdx, 2})) >= 2
            merged_contours{end+1} = unique([contour_pairs{jdx, 1}, contour_pairs{jdx, 2}]);
            merged_idx = unique([merged_idx, idx(jdx, :)]);
            n_merges = n_merges + 1;
        end
    end
    % remove merged contours from input
    contours(merged_idx) = [];
    % repeat merging of merged contours
    if n_merges
        merged_contours = merge_contours(merged_contours);
    end
    % add unmerged contours to output
    for idx = 1:numel(contours)
        merged_contours{end+1} = contours{idx};
    end
end

end
