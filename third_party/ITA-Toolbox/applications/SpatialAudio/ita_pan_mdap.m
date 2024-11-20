function [g] = ita_pan_mdap(LS_pos, VS_pos, varargin)
%ITA_PAN_MDAP Multiple Direction Amplitude Panning
% Description in:
%   [1] Pulkki, Uniform spreading of amplitude panned virtual sources, 1999
% Algorithm creates the spreading sources on a spherical cap around the
% virtual source position. The sources are distributed equally on the cap's
% perimeter which is at an angle of alphaSpread from the caps pole, i.e. 
% the virtual sound source.
% Arguments:
%   LS_pos: itaCoordinates, loudspeaker setup
%   VS_pos: itaCoordinates, virtual source position
% Options: (default)
%   {'spreadAngle'    : (10)  , angular offset between virtual and spread sources
%    'spreadSourceNum': (20)  , number of spreading sources
%    'distanceLoss'   : (true), include 1/r - law
%                       false , panning only}
%   for more Options see ita_pan_vbap.m
%
% Returns:
%   g: MDAP weights

%% Parse input
params.spreadAngle = 10;
params.spreadSourceNum = 20;
params.method = 'static';
params.distanceLoss = true;
[params, ~] = ita_parse_arguments(params, varargin);

%% Calculate spreaded sources
g = zeros(VS_pos.nPoints, LS_pos.nPoints);
VS_spread_base = itaCoordinates(params.spreadSourceNum);
VS_spread = itaCoordinates(params.spreadSourceNum);
% Rotation matrices
X = @(phi) [[1 0 0];[0 cos(phi) -sin(phi)];[0 sin(phi) cos(phi)]];
Y = @(phi) [[cos(phi) 0 -sin(phi)];[0 1 0];[sin(phi) 0 cos(phi)]];
Z = @(phi) [[cos(phi) -sin(phi) 0];[sin(phi) cos(phi) 0];[0 0 1]];
% Precalculate Spread sources
p0 = itaCoordinates(VS_pos.r(1)*[cosd(params.spreadAngle), 0, sind(params.spreadAngle)]);
theta_rot = 0:2*pi/params.spreadSourceNum:2*pi*(1-1/params.spreadSourceNum);
for n = 1:params.spreadSourceNum
    VS_spread_base.cart(n, :) = p0.cart * X(theta_rot(n));
end

for idx = 1:VS_pos.nPoints
    % Rotate to VS_pos
    R = Y(VS_pos.theta(idx) - pi/2)*Z(-VS_pos.phi(idx));
    VS_spread.cart = VS_spread_base.cart*R;

    %% Apply VBAP to all virtual sources
    VS_all = itaCoordinates([VS_pos.cart(idx, :); VS_spread.cart]);
    D_vbap = ita_pan_vbap(LS_pos, VS_all, 'method', params.method, 'normalize', false);
%     activity_factor = sqrt(sum(any(D_vbap ~=0, 2))/LS_pos.nPoints);
%     activity_factor = 1/sqrt(sum(any(D_vbap ~=0, 2)));
    activity_factor = 1;

    %% Normalize weights
    g(idx, :) = sum(D_vbap, 2);
    g(idx, :) = activity_factor * abs(g(idx, :))/norm(g(idx, :));

    if params.distanceLoss
        g(idx, :) = LS_pos.r.'/VS_pos.r(idx) .* g(idx, :);
    end
end

% transpose for easier usage
g = g.';

end
