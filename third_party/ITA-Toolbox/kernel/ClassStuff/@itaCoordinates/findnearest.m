function [ind, dist] = findnearest(this,coords,system,num)

% [ind, dist] = findnearest(this,coords,system,num)
%
%   option: (default)
%       coords () - coordinates for which to find nearest neighbour
%                   can be either itaCoordinate or nx3 Vector
%       system ('cart') - the coordinate system of coords if a vector is
%                         provided 
%       num (1)   - number of neighbours to return


% <ITA-Toolbox>
% This file is part of the ITA-Toolbox. Some rights reserved.
% You can find the license for this m-file in the license.txt file in the ITA-Toolbox folder.
% </ITA-Toolbox>


if ~exist('system','var') || isempty(system)
    % assume cartesian coordinates if no system is given
    system = 'cart';
end

if ~exist('num','var')
    %return one nearest if not otherwise asked
    num = 1;
end

if isnumeric(coords)
    coords = itaCoordinates(coords,system);
end

%speedup for euclidian distance calculations
coordsCart = makeCart(coords);
thisCart = makeCart(this);

if exist('KNNSearch','file') == 3 && ~isempty(thisCart.mPtrtree) % Only if external mex file exists and there are a lot of elements
    %% Using external mex file
    [ind,dist] = KNNSearch(thisCart.cart,coordsCart.cart,thisCart.mPtrtree,num);
else
    %% Old one, using Matlab code
    for idinput = 1:coordsCart.nPoints
        %get euclidian distance to current point
        dists = sqrt(sum((thisCart.cart-repmat(coordsCart.cart(idinput,:),size(thisCart.cart,1),1)).^2,2));
        
        for idx = 1:num
            [dist(idinput, idx), ind(idinput, idx)] = min(dists); %#ok<AGROW>
            
            if sum(dists == min(dists)) > 1 
                % we found multiple potenial points... to pick one, take 
                % angles of potential matches into account
                potentialMatchIdx = find(dists == min(dists));
                
                distsSph = sqrt(sum(this.sph(potentialMatchIdx,2:3)...
                    -repmat(coords.sph(idinput,2:3), numel(potentialMatchIdx),1).^2,2));
                
                %overwrite with better match
                [tmpDist,tmpIdx]   = min(distsSph);
                dist(idinput, idx) = tmpDist;
                ind(idinput, idx)  = potentialMatchIdx(tmpIdx);
            end
            %overwrite distance of current match to exlude in susequent
            %search needed to return n nearest neighbours
            dists(ind(idinput, idx)) = inf;
        end
    end
end
end
