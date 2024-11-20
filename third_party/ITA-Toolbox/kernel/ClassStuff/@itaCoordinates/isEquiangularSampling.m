function [isEquiangular, range] = isEquiangularSampling(this,varargin)
%ISEQUIANGULARSAMPLING(options): check if provided coordinates
%are of an equiangular sampling
%
% check if angular difference of unique theta and phi angles is constant
% and the all points in the meshgrid set by the two are contained in the
% coordinates object
%
%   options (default):
%       tol (1e-6) : angular tolerance in rad when comparing angles
%
% Hark Braren -- hark.braren@akustik.rwth-aachen.de
% 07.12.21

% <ITA-Toolbox>
% This file is part of the ITA-Toolbox. Some rights reserved. 
% You can find the license for this m-file in the license.txt file in the ITA-Toolbox folder. 
% </ITA-Toolbox>

%% initialize and parse input arguments
sArgs.tol = 1e-6;

sArgs = ita_parse_arguments(sArgs,varargin);

tol = sArgs.tol;

isEquiangular  = false; %initialize

%% check data
[thetaU] = unique(this.theta);
meanAngleDelta_theta = mean(diff(sort(thetaU)));

[phiU] = unique(this.phi);
meanAngleDelta_phi = mean(diff(sort(phiU)));


thetaSamplingIsEquiangular = all(diff(sort(thetaU))- meanAngleDelta_theta < tol);
phiSamplingIsEquiangular   = all(diff(sort(phiU))  - meanAngleDelta_phi   < tol);

%% some additional info
range.theta = [min(thetaU),max(thetaU)];
range.phi   = [min(phiU),max(phiU)];

%% checks
if ~thetaSamplingIsEquiangular || ~phiSamplingIsEquiangular
    isEquiangular  = false;
    return;
end

[thetaGrid,phiGrid] = meshgrid(thetaU,phiU);
anglesMissing = setdiff([this.theta,this.phi],[thetaGrid(:),phiGrid(:)],'sorted','rows');

if isempty(anglesMissing)
    isEquiangular = true;
end





