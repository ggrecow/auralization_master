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

%% ---ART Tutorial - Overview---
%  -----------------------------
% This is a basic example for getting started with the ART framework.
% It describes how to setup the Ray Tracing model with an atmosphere.
% Then, in the first part, an example for tracing rays using a defined set
% of launch directions is given. In the second part, the ART framework is
% used to find eigenrays between a source and a receiver. Finally, it is
% shown how the accuracy of eigenrays can be checked.


%% ---Checking ARTMatlab version---
%  --------------------------------
AtmosphericRayTracer.Info();


%% ---Defininition of atmosphere---
%  --------------------------------
atmos = StratifiedAtmosphere;

%Settings
%(These are the default settings, StratifiedAtmosphere is already initialized with these values)
atmos.windProfile = 'log';          %String: 'zero', 'constant', 'log'
atmos.temperatureProfile = 'isa';   %String: 'constant', 'isa'
atmos.humidityProfile = 'constant'; %String: 'constant'

atmos.surfaceRoughness = 0.1;       %Surface Roughness for Log Wind Profile [m]
atmos.frictionVelocity = 0.6;       %Friction velocity for Log Wind Profile [m/s]

atmos.constWindDirection = [1 0 0]; %Normal in wind direction []

atmos.constTemperature = 293;       %Temperature for Const mode [K]
atmos.constStaticPressure = 101325; %Static pressure used in a constant temperature profile

atmos.constWindVelocity = 20;       %Wind velocity for Const Wind Profile [m/s]

atmos.constRelHumidity = 50;        %Constant Realitive Humidity [%]




%% -----Trace Rays-----
%  --------------------
art = AtmosphericRayTracer;

%% Source and initial directions
source = [0 0 100];

deltaAngle = 30;
phi0 = 0:deltaAngle:360-deltaAngle;
theta0 = 90:deltaAngle:180-deltaAngle;
[theta0, phi0] = meshgrid(theta0, phi0);
theta0 = [theta0(:); 180]; %adding south pole
phi0 = [phi0(:); 0];       %adding south pole

rayNormalVectors = AtmosphericRayTracer.AnglesToNormalVectors(theta0, phi0);

%% Ray tracing settings
%--General--
art.maxPropagationTime = 10;                    %Maximum propagation time of rays [s]
art.integrationTimeStep = 0.1;                  %Base time-step for integration of rays [s]
art.solverMethod = 'runge-kutta';               %Method for solving the ODE ('runge-kutta' or 'euler')

%% Run
rays = art.TraceRays(atmos, source, rayNormalVectors);

%% Plot
rays.plot();




%% ---Finding Eigenrays---
%  -----------------------
art = AtmosphericRayTracer;

%% Source and receiver setup
source = [-1000 0 20];
receiver = [0 0 20];

%% Eigenray search settings
%--General--
art.maxPropagationTime = 15;            %Maximum propagation time of rays [s]
art.maxReflectionOrder = 1;             %Maximum number of considered reflections (=number of eigenrays) [] (integer)
art.abortOnReceiverDistIncrease = true; %If enabled, ray tracing will be aborted as soon as ray receiver distance increases (boolean)

%--Adaptive ray zooming--
art.maxReceiverRadius = 1;              %Maximum value for receiver radius [m]
art.maxSourceReceiverAngle = 1;         %Maximum allowed angle between source and receiver sphere [°]
art.maxAngleForGeomSpreading = 0.01;    %Maximimum delta angle of initial direction of neighboring rays used for the calculation of the spreading loss [°]

art.abortMaxNAdaptations = 30;          %Abort criterion for adaptive ray zooming: maximum number of iterations [] (integer)
art.abortMinAngleResolution = 0.001;    %Abort criterion for adaptive ray zooming: minimum angular resolution of rays [°]

%% Run
eigenrays = art.FindEigenrays(atmos, source, receiver);

%% Checking for accuracy of results
% isEigenray property is positive, since rays are result of an eigenray search
eigenrays(1).isEigenray
% The defined receiver sphere was hit and therefore the following property is true
eigenrays(1).receiverSphereHit

%% Plot
eigenrays.plot();
plot3(receiver(1), receiver(2), receiver(3), ' ok', 'MarkerFaceColor', [0 0.6 0]);

%% Performance meta data
% FindEigenrays provides the total number of traced rays as optional second return value
[eigenrays, totalNumRays] = art.FindEigenrays(atmos, source, receiver);
% Each eigenray provides the number of iterations used during its ray zooming process
eigenrays(1).rayZoomingIterations;
eigenrays(2).rayZoomingIterations;

fprintf('Total number of traced rays...\n- Direct path: %i\n- Ground reflection: %i\nNumber of iterations used during adaptive ray zooming process...\n- Direct path: %i\n- Ground reflection: %i\n', totalNumRays(1), totalNumRays(2), eigenrays(1).rayZoomingIterations, eigenrays(2).rayZoomingIterations);



%% ---Limitation of finding eigenrays---
%  -------------------------------------
% In certain cases, finding an eigenray is not possible. A typical example
% is when the receiver resides within the so-called shadow zone which
% cannot be entered by rays. To filter those cases, the receiverSphereHit
% property indicates whether the predefined receiver sphere was hit or not.

%% Source and receiver setup
% The following receiver resides with the shadow zone which is cause by
% strong upward refraction due to upwind conditions
source = [1000 0 20];
receiver = [0 0 1.8];

%% Run
eigenrays = art.FindEigenrays(atmos, source, receiver);

%% Checking for accuracy of results
% isEigenray property is positive, since rays are result of an eigenray search
eigenrays(1).isEigenray
% The defined receiver sphere was NOT hit and therefore the following property is false
eigenrays(1).receiverSphereHit

%% Plot
eigenrays.plot();
plot3(receiver(1), receiver(2), receiver(3), ' ok', 'MarkerFaceColor', [0 0.6 0]);
view(0,0);