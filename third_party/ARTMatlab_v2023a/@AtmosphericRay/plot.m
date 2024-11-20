function varargout = plot(rays, varargin)
%Plots ray(s) and given source(s).
%
%Syntax is similar to Matlab plot functions:
%If the first of the inputs is an axes handle, everything is plotted into
%these axes. If not, a new figure is opened and axes are created. Besides
%that, the function expects the typical sequence of property names and values.
%
%Inputs:
%Modes:
%Property names (default values):
%'LaunchElevation' (none):  If passing this string, ray color corresponds
%                           to launch elevation angle
%'SourceColor' ([0 0 1]):   RGB vector with source color.
%                           If empty, source(s) won't be plotted.
%Line plot properties:      Additional inputs are fed directly to plot3 function
%                           (type 'doc line properties' for more information).
%                       
%Output:
%ax:                    Axes handle containing all plots
%
%Examples:
%   rays    - array of AtmosphericRay
%   ax      - axes handle
%   rays.plot(ax, 'SourceColor', [1 0 0], 'LineWidth', 2, 'LineStyle', '--');
%   
%   rays.plot('Color', [0 1 0]);
%   
%   rays.plot(ax, 'LaunchElevation', 'SourceColor', [1 0 0], 'LineWidth', 2);

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


%% Input Checks
assert(nargout <= 1, 'Only one output argument allowed')
assert(~isempty(rays) && rays(1).numPoints > 0, 'AtmosphericRay is empty. Nothing to plot.')

sourceColor = [0 0 1];
idxSourceColor = find( strcmpi(varargin, 'SourceColor') );
assert(numel(idxSourceColor) <= 1, 'Found duplicated option "SourceColor"')

if numel(idxSourceColor) == 1
    assert(numel(varargin) >= idxSourceColor+1, 'No value for "SourceColor" property found')
    sourceColor = varargin{idxSourceColor+1};
    varargin([idxSourceColor idxSourceColor+1]) = [];
end
assert(isempty(sourceColor) || isnumeric(sourceColor) && numel(sourceColor) == 3, 'Option SourceColor must be an rgb vector or empty')

%% Axes given?
if numel(varargin) && isa(varargin{1}, 'matlab.graphics.axis.Axes')
    ax = varargin{1};
    varargin = varargin(2:end);
else
    ax = axes(figure);
end
hold(ax, 'on')

%% Launch Elevation Mode?
idx = cellfun(@(x)ischar(x)&&strcmp(x,'LaunchElevation'),varargin);
bLaunchElevation = any(idx);
varargin(idx) = [];

%% Plots
sources = zeros(numel(rays), 3);
for idx=1:numel(rays)
    plotSingleRay(ax, rays(idx), bLaunchElevation, varargin{:});
    sources(idx,:) = rays(idx).r0.cart;
end

if ~isempty(sourceColor)
    uniqueSources = unique(sources,'rows');
    for idx=1:size(uniqueSources,1)
        plotSource(ax, uniqueSources(idx,:), sourceColor);
    end
end

%% Plot Layout
xlabel(ax, 'x [m]')
ylabel(ax, 'y [m]')
zlabel(ax, 'z [m]')
view(ax, 30, 30)
grid(ax, 'on')
rotate3d(ax, 'on')

%% Lauch Elevation mode plot layout
if bLaunchElevation
    [colorMap, thetaDeg] = LaunchElevationColorMap();
    colormap(ax, colorMap);
    ax.CLim = [thetaDeg(1) thetaDeg(end)];
    c = colorbar(ax);
    c.Label.String = 'Launch elevation [°]';
    c.Ticks = -90:30:90;
end

%% Output
if nargout; varargout{1} = ax; end


function plotSingleRay(ax, ray, bLaunchElevation, varargin)

if(ray.numPoints < 1)
    error('AtmosphericRay: Cannot plot empty ray.')
end

[r, ~] = mirrorGroundReflection(ray);

%% LaunchElevationMode
if bLaunchElevation
    color = getLaunchElevationColor(ray);
    varargin = [varargin, {'Color', color}];
end

%% Ray Path
plot3(ax, r.x, r.y, r.z, varargin{:})

function plotSource(ax, source, sourceColor)
plot3(ax, source(1), source(2), source(3), ' ok', 'markerfacecolor', sourceColor);

function [r, n] = mirrorGroundReflection(ray)

%Mirrors all points and wavefront normals which are below the ground at the
%plane z=0.

r = ray.r;
n = ray.n;

idxBelowGround = r.z<0;

r.z = abs(r.z);
n.z(idxBelowGround) = -n.z(idxBelowGround);


function color = getLaunchElevationColor(ray)
[colorMap, thetaDeg] = LaunchElevationColorMap();

launchElevationDeg = 90-ray.n0.theta_deg;
color = interp1(thetaDeg, colorMap, launchElevationDeg);

function [colorMap, thetaDeg] = LaunchElevationColorMap()
colorMap = parula();

thetaDeg = linspace(-90, 90, size(colorMap, 1));
