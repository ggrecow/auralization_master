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

%% ---ART Tutorial - Propagation Model---
%  --------------------------------------
% The Atmospheric Ray Tracing framework allows to calculate acoustic
% propagation parameters based on simulated eigenrays. Those parameters can
% be used for an auralization process (auralization parameters).
% Additionally, it is possible to calculate a transfer function describing
% the sound transmission for a static source and receiver.
%
% This example introduces the AtmosphericPropagation class which is
% designed for those purposes.


%% ---Calculate Eigenrays---
%  -------------------------
atmos = StratifiedAtmosphere;
art = AtmosphericRayTracer;
source = [-1000 0 20];
receiver = [0 0 20];
eigenrays = art.FindEigenrays(atmos, source, receiver);


%% ---Initializing Propagation Model---
%  ------------------------------------
propagationModel = AtmosphericPropagation;

%% Atmosphere
% It is important to use the same atmosphere as for the eigenray
% simulation.
propagationModel.atmosphere = atmos;

%% Frequency vector
samplingRate = 44100;
tMax = 3;
nBins = samplingRate*tMax/2;
propagationModel.frequencyVector = linspace(0, samplingRate/2, nBins);


%% Setting the reflection factor
% The ground reflection factor is required for the transfer function
% calculation. It is 1 by default. It can either be a singe value or a
% complex-valued vector with same length as the frequency vector.
propagationModel.groundReflectionFactor = 0.9 * ones(size(propagationModel.frequencyVector));



%% ---Auralization parameters---
%  -----------------------------

%% Geometrical spreading
% The spreading loss factor is calculated at the end of the eigenray
% search based on the Blokhintzev invariant. It is stored in a property
% of the eigenray.
spreadingLossDirect = eigenrays(1).spreadingLoss
spreadingLossReflection = eigenrays(2).spreadingLoss

%% Propagation delay
% The propagation delay is implicitely stored in the time vector of the
% eigenrays
tDelayDirect = eigenrays(1).t(end)      %[s]
tDelayReflection = eigenrays(2).t(end)  %[s]

%% Air attenuation
% The air attenuation is calculated by integrating the attenuation defined
% in ISO 9613-1 [dB/m] along the ray path.
attenuation = propagationModel.AirAttenuation(eigenrays);

%Plot
attenuationITA = itaResult(attenuation, propagationModel.frequencyVector, 'freq');
attenuationITA.comment = 'Air attenuation';
attenuationITA.channelNames = {'Direct path', 'Reflection'};
attenuationITA.pf();



%% ---Transfer function---
%  -----------------------
% The transfer function (TF) is calculated for each eigenray individually and
% then combined using the priciple of superposition.
% It considers spreading loss, propagation delay, air attenuation and a
% reflection factor for the ground reflection.

%% Calculation
[combinedTF, individualTFs] = propagationModel.TransferFunction(eigenrays);

%% Plot
allTFs = ita_merge(individualTFs, combinedTF);
allTFs.channelNames = {'Direct path TF', 'Ground reflection TF', 'Combined TF'};
allTFs.pf();
