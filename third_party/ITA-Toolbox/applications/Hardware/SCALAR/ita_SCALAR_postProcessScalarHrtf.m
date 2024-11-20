function resultHRTF = ita_SCALAR_postProcessScalarHrtf(HRTF_raw,options)
%ITA_SCALAR_POSTPROCESSSCALARHRTF - +++ Short Description here +++
%  This function ++++ FILL IN INFO HERE +++
%
%  Syntax:
%   audioObjOut = ita_SCALAR_postProcessScalarHrtf(audioObjIn, options)
%
%   options struct is mandatory containing:
%       options.ref = 'pathToRef.ita' path to reference measurement
%       options.tw  = [start stop] ita_time_window option
%
%  Example:
%   audioObjOut = ita_SCALAR_postProcessScalarHrtf(audioObjIn)
%
%  See also:
%   ita_SCALAR_measureHRTF,
%
%   Reference page in Help browser 
%        <a href="matlab:doc ita_SCALAR_postProcessScalarHrtf">doc ita_SCALAR_postProcessScalarHrtf</a>

% <ITA-Toolbox>
% This file is part of the ITA-Toolbox. Some rights reserved. 
% You can find the license for this m-file in the license.txt file in the ITA-Toolbox folder. 
% </ITA-Toolbox>


% Author: Hark Braren -- Email: hark.braren@akustik.rwth-aachen.de
% Created:  21-Sep-2022 


%% get reference
ref = ita_read(options.ref);
% ref = merge(ref);
ref_tw =   ita_time_window(ref,options.tw,'time');
% make end sample be div by 4 for daff export
endSample = round(options.tw(2) .* ref(1).samplingRate)+1;
endSample = endSample + mod(endSample,4);
ref_tw      =    ita_time_crop(ref_tw,[1 endSample],'samples');


%#ok<*AGROW> 
for iCh = 1:2
    ref_finished(iCh) = ita_smooth_notches(ref_tw(iCh),'bandwidth',1/2,...
        'threshold', 3); 
    ref_finished(iCh) = ita_time_shift(ref_finished(iCh),-2e-3,'time');
end


%% Process the measurement
%1) time window and cropping
result_tw   =  ita_time_window(HRTF_raw,options.tw,'time');
result_crop =  ita_time_crop(result_tw,[1 endSample],'samples');

% divide by reference
for iCh = 1:2
    data_referenced(iCh) = ita_divide_spk(result_crop(iCh),ref_finished(iCh),'regularization',[100 20000]);
end


%% smooth low frequencies
allMeasurements_full = ita_HRTF_postProcessing_smoothLowFreq(data_referenced,'cutOffFrequency',200,'upperFrequency',1000,'timeShift',0);


%% add Metadata
for index = 1:length(allMeasurements_full)
    allMeasurements_full(index).userData = options;
end

%% append commit id to history
commitID = ita_git_commit_id;
for index = 1:length(allMeasurements_full)
    if ~isempty(commitID)
        allMeasurements_full(index) = ita_metainfo_add_historyline(allMeasurements_full(index),'ita_SCALAR_postProcessHrtf',commitID);
    end
end

result = allMeasurements_full;

%turn into itaHRTF
resultHRTF = itaHRTF(result);
end