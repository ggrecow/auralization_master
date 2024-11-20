%ITA_SCALAR_MEASUREHRTF - measurement script for an HRTF using the SCALAR array
%  

% <ITA-Toolbox>
% This file is part of the ITA-Toolbox. Some rights reserved. 
% You can find the license for this m-file in the license.txt file in the ITA-Toolbox folder. 
% </ITA-Toolbox>


% Author: Hark Braren -- Email: hark.braren@akustik.rwth-aachen.de
% Created:  25-Jun-2021 


%set to true if you want to confirm the latency measurement. PLease be sure what you are doing !!!
options.measureLatency = false;
%get default loudspeaker positions !! BE AWARE OF UNDOCUMENTED CHANGES AND MOVED LOUDSPEAKERS!!
speakerLocations = ita_SCALAR_loudspeakerPositions('type','ideal');

subjID = 999;

saveFolder = fullfile('./'); %Where to save the HRTF
saveName   = sprintf('Participant_%d_%s',subjID,datestr(now,'dd-mmm-yy'));

if ~exist(saveFolder,'dir')
   mkdir(saveFolder) 
end

%% measurement Setup
ms = itaMSTFinterleaved;

if options.measureLatency == false
    ms.latencysamples = 1724;
else
    %% measure Latency
    ms.inputChannels = 2;
    ms.outputChannels = 43;
    
    input('Connect the amp to ch 2')
    ms.run_latency;
    input('Disconnect the amp')
end


ms.outputChannels = ([1:68,79:84]); %update to new channels
ms.inputChannels  = [5,6];
ms.freqRange = [100 22050];
ms.outputamplification = -20;

% prepare for measurement
ms.optimize('mode','standard','plot',0,'sweeprate_range',[2 6])

%% run the measurement

result = ms.run;

% prepare for post processign and save raw data
result = [result.ch(1).merge,result.ch(2).merge];
result(1).channelCoordinates = speakerLocations;
result(2).channelCoordinates = speakerLocations;
ita_write_ita(result,fullfile(saveFolder,[saveName,'_raw.ita']));


%% post Processing - ch
% refMeasurement = ''; %path to the reference Measurement
% options.ref = refMeasurement;
% options.tw = [0.0055, 0.006];
% result_pp = ita_SCALAR_postProcessScalarHRTF(result,options);
% ita_write_ita(result_pp,fullfile(saveFolder,[saveName,'.ita']));