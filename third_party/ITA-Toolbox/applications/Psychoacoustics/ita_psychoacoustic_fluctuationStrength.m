function varargout = ita_psychoacoustic_fluctuationStrength(varargin)
%ITA_CALCULATE_FLUCTUATIONSTRENGTH - +++ Short Description here +++
%  This function ++++ FILL IN INFO HERE +++
%
%  Syntax:
%   [double, itaResultObj, double] = ita_calculate_fluctuationStrength(audioObjIn, options)
%
%   Options (default):
%           'modFreq' ('auto-detect')        : modulation frequency (if known: scalar or two-element vector between [0.1,100])
%           'soundFieldType' ('free')        : field conditions ('diffuse')
%
%  Example:
%   [fluc_tot, fluc_spec, f_mod] = ita_calculate_fluctuationStrength(audioObjIn, 'modFreq', 20)
%
%  See also:
%   ita_calculate_loudness, ita_calculate_sharpness
%
%   Reference page in Help browser 
%        <a href="matlab:doc ita_calculate_fluctuationStrength">doc ita_calculate_fluctuationStrength</a>

% <ITA-Toolbox>
% This file is part of the ITA-Toolbox. Some rights reserved. 
% You can find the license for this m-file in the license.txt file in the ITA-Toolbox folder. 
% </ITA-Toolbox>


% Author: ShaimaaDoma -- Email: sdo@akustik.rwth-aachen.de
% Created:  05-Oct-2021 


%% Initialization and Input Parsing
sArgs   = struct('pos1_data',       'itaAudio', ...
                 'modFreq',         'auto-detect',... % 'auto-detect', scalar or 2-elem-vector if pre-known
                 'soundFieldType',  'free');          % 'free', 'diffuse'
                                            
[input,sArgs] = ita_parse_arguments(sArgs,varargin); 

%% Checks

% Matlab version
matlabVer = version('-release');
if str2double(matlabVer(1:4))< 2020
    error('ITA_CALCULATE_FLUCTUATIONSTRENGTH: A newer version of Matlab (at least 2020a) is required for running this function. Use ita_loudness alternatively.')
end

%% Core
data_in = input.timeData;  
fs = input.samplingRate;

if strcmp(input.channelUnits, 'Pa')
else
 ita_verbose_info('The expected unit of the itaValue object is Pascal. This does not seem to be the case for your object. Double check because Pascal is assumed as unit during the calculation.', 0)
end 

if nargout == 0
    if size(data_in,2)>2
        for idxAudioCh = 1:size(data_in,2)
            figure
            acousticFluctuation(data_in(:,idxAudioCh),fs,1,'SoundField',sArgs.soundFieldType,...
                'ModulationFrequency', sArgs.modFreq);
        end
    else
        acousticFluctuation(data_in,fs,1,'SoundField',sArgs.soundFieldType,...
            'ModulationFrequency', sArgs.modFreq);
    end

else
    if size(data_in,2)<=2
        % use input data directly
        [fluc, flucS, f_mod] = acousticFluctuation(data_in,fs,1,'SoundField',sArgs.soundFieldType,...
            'ModulationFrequency', sArgs.modFreq);
    else
        fluc = []; flucS = [];
        for idxAudioCh = 1:size(data_in,2)
            [fluc_temp, flucS_temp, f_mod_temp] = acousticFluctuation(data_in(:,idxAudioCh),fs,1,'SoundField',sArgs.soundFieldType,...
                    'ModulationFrequency', sArgs.modFreq);
            fluc(:,idxAudioCh) = fluc_temp;
            flucS(:,:,idxAudioCh) = flucS_temp;
            f_mod(idxAudioCh) = f_mod_temp;
        end
    end
end

% TODO: check why time resolution is not accepted yet as an argument.

%% Set Output
if nargout > 0
    totalFluc = fluc;
    varargout(1) = {mean(totalFluc(500:end,:))}; %first 500 frames are excluded (see MATLAB Doku)
    
    if nargout > 1
        freqVector         = (0.5: 0.5: 23.5)';
        timeRes            = 2e-3; % in acousticFluctuation only 2ms as time resolution
        specFluc                    = itaSpecificResult(flucS, 0:timeRes:(numel(flucS(:,1))-1)*timeRes, freqVector); % copy audio struct information from input struct
        for idxAudioCh = 1:size(data_in,2)
            specFluc.channelUnits(idxAudioCh,:)    = {'vacil/Bark'}; 
            specFluc.channelDataType(idxAudioCh,:)    = {'Specific Fluctuation Strength'};
        end
        varargout(2) = {specFluc};

        if nargout > 2
            varargout(3) = {f_mod};
        end
    end

end