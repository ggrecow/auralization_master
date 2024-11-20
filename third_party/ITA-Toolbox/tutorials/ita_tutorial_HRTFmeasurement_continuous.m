%% Continuous HRTF measurement using the HRTFarc 
%
% This tutorial demonstrates the use of the ITA Arc for a continuous HRTF
% measurement, where the participant is static and the arc is in motion. It
% is based on the script "ita_HRTFarc_measurementScript_continuous.m" with
% further notes and extended by the reference measurement.
% 
% For information about the presented methods, refer to
%   Richter & Fels - "On the influence of continuous subject rotation during
%   high-resolution HRTF measurements", 2019

% Author: Shaimaa Doma -- Email: sdo@akustik.rwth-aachen.de
% Created:  Oct-2021

% <ITA-Toolbox>
% This file is part of the ITA-Toolbox. All rights reserved.
% You can find the license for this m-file in the application folder.
% </ITA-Toolbox>

%% General settings

fs = 48000;             % TODO: set wanted sampling rate  
ita_preferences('samplingRate',fs)

% TODO in code or manually in ita_preferences pop-up window: set Movtec COM Port to 7  

%% === Prepare reference measurement ===
% Two Sennheiser KE3 microphones are to be measured simultaneously at the
% center of the arc. Each of the 64 loudspeakers emits a single sweep, all
% together forming an interleaved signal. The arc is stationary.

% init measurement object
ms = itaMSTFinterleaved; 
ms.latencysamples = 629;        % used per default, unless another latency measurement 
                                % is performed via ms.run_latency
                                
% I/O
ms.outputChannels = (1:64)+18;  % 64 loudspeakers as output
ms.inputChannels = [12,13];     % TODO: insert correct numbers
                                % (RME Octamic channel numbers +10) for L and R ear microphone
micNumbers = [3,4];             % TODO: insert microphone numbers (according to sticker)

% signal properties
ms.freqRange = [300 fs/2];      % frequency range of the sweeps
ms.samplingRate = ita_preferences('samplingRate');
ms.optimize                     % optimize sweep rate
ms.twait = 0.03;                % delay of 30 ms between subsequent sweeps

                                   
%% === Run reference measurement & check the SNR ===

res = ms.run;

% Take a look
res.merge.pf
res.merge.ptd
% You may need to adjust the gain manually to improve the SNR. 
% In that case, re-run this section !!

ita_write(res,sprintf('ref_octaCh_%d_%d_micCh_%d_%d',...
                        ms.inputChannels(1)-10,ms.inputChannels(2)-10,micNumbers(1),micNumbers(2)))

%% === Prepare for HRTF measurement ===
% The microphones are placed at the ear canal entrances of the participant,
% who is now positioned at the center. The arc rotates around them while
% emitting a continuous signal of interleaved sweeps. 

clear ms

% === Create a measurement object with almost similar settings ===
ms = itaMSTFinterleaved; 
ms.latencysamples = 629;  

% I/O
ms.outputChannels = (1:64)+18; 

% use the same microphones and Octamic channels as before !!
ms.inputChannels = [12,13,14];  % (!!!) include a third channel that receives the mechanical switch signal 
                                % indicating the start and end of a full rotation: [L, R, switch]
                                % The order of the channels is flexible but later needs to be 
                                % considered in post-processing.                              
micNumbers = [3,4];                      

% signal properties
ms.freqRange = [300 fs/2];      
ms.samplingRate = ita_preferences('samplingRate');
ms.optimize                     
ms.twait = 0.03;                

% === Set the arc rotation speed (-->  azimuth resolution) ===
% The number of repetitions refers to the number of sweeps emitted by each
% lourspeaker. The higher the rep number and the resolution, the slower the arc
% movement and the longer the measurement duration (64 ~ 3 minutes).
numRepetitions = 144;           % e.g. 360°/2.5° = 144 corresponds to an azimuth resolution of 2.5°. 
ms.repetitions = numRepetitions;

% === Prepare a motor object ===

% init
iMS = test_itaEimarMotorControl;
% The connected motor 

% hand the measurement object to motor object
iMS.measurementSetup = ms;

% back up (exact excitation signals might be needed for manual deconvolution in case of unexpected issues)
save('setupSave.mat','iMS','ms');

%% === Run the HRTF measurement ===

% make reference move, then rotate back by 45° for the starting position
iMS.prepareContinuousMeasurement;                                      

% actual measurement
[res,raw] = iMS.runContinuousMeasurement;

ita_write(res,'measConti')      % after deconvolution with sweep
ita_write(raw,'measConti_raw')  % before deconvolution 

% Take a look. You should see the head shadowing variations for both
% ears, as well as switch signal peaks close to the beginning and end of the
% measurement period.
res_raw.pt

% move the arc back (faster than during actual measurement)
iMS.moveTo('HRTFArc',20, 'absolut', true, 'speed', 5, 'wait', true);   

% ALWAYS leave the arc in reference position !!!
iMS.reference

%% Clean up
% clear workspace (or at least the measurement and motor objects) before
% attempting attempting another run

% ccx
% clear ms iMS

%% In case of emergency
% iMS.stop



