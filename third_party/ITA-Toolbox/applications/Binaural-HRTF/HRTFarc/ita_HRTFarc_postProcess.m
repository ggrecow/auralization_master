function result = ita_HRTFarc_postProcess(varargin)
%ITA_HRTF_ARC_POSTPROCESS - postprocessing of stepwise HRTF measurements
%  This function is used to process HRTF measurements created using
%  ita_HRTFarc_measurementScript. It takes care of time windowing the
%  measurement and reference files, and spectral division by the reference
%  meaurement. Additional options such as centering of azimuth angle based
%  on ITD are available (see options and code)
%  The returned object is an 1x2 itaAudio with the HRTF measurements for
%  each ear. To get an HRTF use: itaHRTF(result)
%
%  Syntax:
%   audioObjOut = ita_HRTF_arc_postProcess(varargin)
%
%   Options (default):
%           'dataFolder'    ('')    : location of the data folder from the measurement 
%           'refFile'       ('')    : location of an ita audio file with reference 
%           'tw'  ([0.006 0.008])   : the time window edges - see ita_time_window  
%           'shiftSamples' (128)    : final hrtfs are circularly shifted by this number of samples to make them causal
%           'dataChannel'   (1:2)   :microphone channels
%           'rotationDirection'  (-1) : rotation direction
%           'normalize'     ([])      : normalization factors
%           'ms'            ([])      : needed only in continous measurements 
%           'adjustAzimzimuth' (false): adjust azimuthAngle based on ITD
%           'itdMethod'     ('xcorr') : method for itd calculation
%           'regularization'([100 20000])) : frequency parameter for regularized inverse
%
%
%  Example:
%   audioObjOut = ita_HRTFarc_postProcess('dataFolder',yourDataFolder,'refFile',locationOfYourRefFile)
%
%  See also:
%   ita_toolbox_gui, ita_read, ita_write, ita_generate
%
%   Reference page in Help browser
%        <a href="matlab:doc ita_HRTF_arc_postProcess">doc ita_HRTF_arc_postProcess</a>

% <ITA-Toolbox>
% This file is part of the ITA-Toolbox. Some rights reserved.
% You can find the license for this m-file in the license.txt file in the ITA-Toolbox folder.
% </ITA-Toolbox>


% Author: Jan-Gerrit Richter -- Email: jri@akustik.rwth-aachen.de
% Created:  05-Oct-2018
% Updated:  25-Jul-2022 (HBR)


%% Initialization and Input Parsing
% all fixed inputs get fieldnames with posX_* and dataTyp
% optonal inputs get a default value ('comment','test', 'opt1', true)
% please see the documentation for more details
sArgs  = struct(...
    'dataFolder','',  ...
    'refFile','',...
    'ms',[], ...
    'tw',[0.006 0.008], ...
    'shiftSamples',[128],...
    'rotationDirection',-1, ...
    'normalize',[], ...
    'adjustAzimzimuth',false,...
    'itdMethod','xcorr', ...
    'dataChannel',1:2,...
    'regularization',[100 20000]);

[options] = ita_parse_arguments(sArgs,varargin);

%% first, load and prepare the reference
ref_load = ita_read(options.refFile);
ref = merge(ref_load);

ref_tw      =   ita_time_window(ref,options.tw,'time');
% make end sample be div by 4 for daff export
endSample = round(options.tw(2) .* ref(1).samplingRate)+1;
endSample = endSample + mod(endSample,4);
ref_tw      =    ita_time_crop(ref_tw,[1 endSample],'samples');



ref_finished = ita_smooth_notches(ref_tw,'bandwidth',1/2,...
    'threshold', 3);

clear tmp;
% create a multi instance again - create
% 1xnMicrophonesInReferenceMEasurement itaAudio with nLS channels
if ref_load(1).nChannels > 1
    for iMicrophone = 1:ref_load(1).nChannels
        ref_finished_perMicrophone(iMicrophone) = ref_finished.ch(iMicrophone:ref_load(1).nChannels:ref_finished.nChannels);
    end
    ref_finished = ref_finished_perMicrophone;
end

%% load all measurements
files = dir(fullfile(options.dataFolder, 'data', '*.ita')); %get all .ita Files
idealCoords = ita_HRTFarc_returnIdealNewArcCoordinates;
wb = itaWaitbar(numel(files));
%load all files
for iFileToLoad = 1:numel(files)
    data = ita_read_ita(fullfile(options.dataFolder, 'data', [num2str(iFileToLoad) '.ita']));
    if length(data) == 1
        data_tmp = options.ms.crop(data);
    else
        data_tmp = data;
    end
    
    for dataIndex = 1:length(options.dataChannel)
        data(dataIndex) = merge(data_tmp.ch(dataIndex));
    end
    
    for iMicrophone = 1:length(options.dataChannel) %go throuh each inear microphone
        tmp = data(iMicrophone);
        coordinates = tmp.channelCoordinates;
        if options.rotationDirection == -1
            phi = mod(0 - mean(unique(coordinates.phi)),2*pi);
        else
            phi = mod(mean(unique(coordinates.phi)),2*pi);
        end
        coordinates = idealCoords;
        coordinates.phi = phi;
        
        tmp.channelCoordinates = coordinates;
        allMeasurementsRaw(iFileToLoad,iMicrophone) = tmp;
        
        
        % Fancy Cropping as for reference
        data_crop      =   ita_time_window(tmp,options.tw,'time');% Q: options.tw = [6ms 8ms]
        % make end sample be div by 4 for daff export
        endSample = round(options.tw(2) .* tmp(1).samplingRate)+1;
        endSample = endSample + mod(endSample,4);
        data_crop      =    ita_time_crop(data_crop,[1 endSample],'samples');
        
        % divide by reference
        allMeasurements(iFileToLoad,iMicrophone) = ita_divide_spk(data_crop,ref_finished(iMicrophone),'regularization',options.regularization);
    end
    wb.inc;
end
for dataIndex = 1:length(options.dataChannel)
    allMeasurements_full(dataIndex) = merge(allMeasurements(:,dataIndex));
end

wb.close;

%% additional postprocessing
% time shift, itd correction and normalization
% shift analog zu ramona
for indexRefFinished = 1:2
    allMeasurements_full(indexRefFinished) = ita_time_shift(allMeasurements_full(indexRefFinished) , sArgs.shiftSamples, 'samples');
end


tmpCoords = allMeasurements_full(1).channelCoordinates;
% calculate ITD and shift to 0 -- search for "ITD == 0"
if options.adjustAzimzimuth
    [centerPoint,itdData] = ita_HRTFarc_pp_itdInterpolate(allMeasurements_full,tmpCoords,options);
    if itdData.error > 0.01
        disp('warning: itd match does not look good. something is wrong in either the data, or the itd method');
    end
    
    for dataIndex = 1:length(options.dataChannel)
        tmpCoords = allMeasurements_full(dataIndex).channelCoordinates;
        tmpCoords.phi_deg = mod(tmpCoords.phi_deg -centerPoint,360);
        allMeasurements_full(dataIndex).channelCoordinates = tmpCoords;
    end
end


% normalize
if ~isempty(options.normalize)
    allMeasurements_full(1).freqData = allMeasurements_full(1).freqData ./ options.normalize(1);
    allMeasurements_full(2).freqData = allMeasurements_full(2).freqData ./ options.normalize(2);
end

allMeasurements_full = ita_HRTF_postProcessing_smoothLowFreq(allMeasurements_full,'cutOffFrequency',500,'upperFrequency',1200,'timeShift',0);

%% add metadata to result
% append options to userdata
for indexRefFinished = 1:length(allMeasurements_full)
    allMeasurements_full(indexRefFinished).userData = options;
end

% append commit id to history
commitID = ita_git_commit_id();
for indexRefFinished = 1:length(allMeasurements_full)
    if ~isempty(commitID)
        allMeasurements_full(indexRefFinished) = ita_metainfo_add_historyline(allMeasurements_full(indexRefFinished),'ita_HRTFarc_postProcessContinuous',commitID);
    end
end

result = allMeasurements_full;
end