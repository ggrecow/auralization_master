function varargout = ita_calculate_roughness(varargin)
%ITA_CALCULATE_ROUGHNESS
%  This is a wrapper function for acousticRoughness.m, which determines
%  total roughness (R) and specific loudness (R') values for one or two
%  audio channels. It applies 'ISO 532-1' after Zwicker for mono or stereo
%  calculations.
%  The output is time-dependent, providing values per 0.5 ms interval.
%  A minimum track duration of 0.5 seconds is recommended. If not priorly
%  known, the most dominant modulation frequency is estimated for the
%  subsequent calculations.
%  Requirement: Matlab 2020a or later.
%
%  Syntax:
%   [double, itaResultObj, double] = ita_calculate_roughness(audioObjIn, options)
%
%   Options (default):
%           'modFreq' ('auto-detect')        : modulation frequency (if known: scalar or two-element vector between [1,1000])
%           'soundFieldType' ('free')        : field conditions ('diffuse')
%
%  Example:
%   [R_tot, R_spec, f_mod] = ita_calculate_roughness(audioObjIn, 'soundFieldType','diffuse')
%
%  See also:
%   ita_calculate_loudness, ita_calculate_sharpness
%
%   Reference page in Help browser
%        <a href="matlab:doc ita_calculate_roughness">doc ita_calculate_roughness</a>

% <ITA-Toolbox>
% This file is part of the ITA-Toolbox. Some rights reserved.
% You can find the license for this m-file in the license.txt file in the ITA-Toolbox folder.
% </ITA-Toolbox>


% Author: ShaimaaDoma -- Email: sdo@akustik.rwth-aachen.de
% Created:  05-Oct-2021

%% Initialization and Input Parsing
sArgs   = struct('pos1_data',       'itaAudio', ...
    'modFreq',         'auto-detect',... % 'auto-detect', scalar or 2-elem-vector if pre-known
    'soundFieldType',  'free', ...        % 'free', 'diffuse'
    'timeResolution','standard');         % 'standard' (2ms) or 'high' (0.5ms)

[input,sArgs] = ita_parse_arguments(sArgs,varargin);

%% Checks

% Matlab version
matlabVer = version('-release');
if str2double(matlabVer(1:4))< 2020
    error('ITA_CALCULATE_Roughness: A newer version of Matlab (at least 2020a) is required for running this function. Use ita_loudness alternatively.')
end

% Check for input unit
if strcmp(input.channelUnits, 'Pa')
else
    ita_verbose_info('The expected unit of the itaValue object is Pascal. This does not seem to be the case for your object. Double check because Pascal is assumed as unit during the calculation.', 0)
end

%% Core
data_in = input.timeData;  % TODO: handling multi-channel (>2) data if needed. Currently only one or two per default.
fs = input.samplingRate;

if nargout == 0
    if size(data_in,2)>2
        for idxAudioCh = 1:size(data_in,2)
            figure
            acousticRoughness(data_in(:,idxAudioCh),fs,1,'SoundField',sArgs.soundFieldType,...
                'ModulationFrequency', sArgs.modFreq);
        end
    else
        acousticRoughness(data_in,fs,1,'SoundField',sArgs.soundFieldType,...
            'ModulationFrequency', sArgs.modFreq);
    end
else
    if size(data_in,2)>2
        R = []; RS = [];
        for idxAudioCh = 1:size(data_in,2)
            [R_temp, RS_temp, f_mod] = acousticRoughness(data_in(:,idxAudioCh),fs,1,'SoundField',sArgs.soundFieldType,...
                'ModulationFrequency', sArgs.modFreq);
            R(:,idxAudioCh) = R_temp;
            RS(:,:,idxAudioCh) = RS_temp;
        end
    else
        [R, RS, f_mod] = acousticRoughness(data_in,fs,1,'SoundField',sArgs.soundFieldType,...
            'ModulationFrequency', sArgs.modFreq);
    end
end

%% Set Output
if nargout > 0

    % totalRoughness = itaValue(R,'aspers');
    totalRoughness = R;
    varargout(1) = {mean(totalRoughness(2001:end, :))}; %only take values from sample 2000 on (see acousticRoughness MATLAB documentation)

    if nargout > 1
        freqVector         = (0.5: 0.5: 23.5)';
        if strcmp(sArgs.timeResolution, 'standard')
            time_res = 2e-3;
        elseif strcmp(sArgs.timeResolution, 'high')
            time_res = 5e-4;
        else
            ita_verbose_info('Something might have changed in the time resolution parameters of the Matlab function acousticLoudness, please check!')
        end
        specRoughness                    = itaSpecificResult(RS, 0:time_res:(numel(RS(:,1))-1)*time_res, freqVector); % copy audio struct information from input struct
        specRoughness.channelUnits(:)    = {'aspers/Bark'};
        for idxAudioCh = 1:size(data_in,2)
            specRoughness.channelUnits(idxAudioCh,:)    = {'aspers/Bark'};
            specRoughness.channelDataType(idxAudioCh,:)    = {'Specific Roughness'};
        end
        varargout(2) = {specRoughness};

        if nargout > 2
            varargout(3) = {f_mod};
        end
    end
end
end