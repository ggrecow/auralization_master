function varargout = ita_calculate_sharpness(varargin)
% ITA_CALCULATE_SHARPNESS - Calculates sharpness from multiple audio channels
%  This is a wrapper function for acousticSharpness.m, which determines 
%  sharpness values for one or two audio channels. Here, more than two
%  channels are admitted.
%  Per default, the function applies DIN 45692 (ISO 532-1) for frequency weighting, 
%  or alternatively 'Aures' or 'von Bismarck' weighting. 
%  Requirement: Matlab 2020a or later.
%
%  Syntax:
%   double = ita_calculate_sharpness(audioObjIn, options)
%
%   Options (default):
%           'weighting' ('DIN 45692')        : frequency weighting('Aures', 'von Bismarck')
%           'soundFieldType' ('free')        : field conditions ('diffuse')
%           'timeVarying' (false)            : time-dependency (true)
%
%  Example:
%   S = ita_calculate_sharpness(audioObjIn, 'weighting', 'von Bismarck')
%
%  See also:
%   ita_toolbox_gui, ita_read, ita_write, ita_generate
%
%   Reference page in Help browser 
%        <a href="matlab:doc ita_calculate_sharpness">doc ita_calculate_sharpness</a>

% <ITA-Toolbox>
% This file is part of the ITA-Toolbox. Some rights reserved. 
% You can find the license for this m-file in the license.txt file in the ITA-Toolbox folder. 
% </ITA-Toolbox>


% Author: ShaimaaDoma -- Email: sdo@akustik.rwth-aachen.de
% Created:  05-Oct-2021 


%% Initialization and Input Parsing
sArgs   = struct('pos1_data',       'itaAudio', ...
                 'weighting',       'DIN 45692',...   % 'DIN 45692', 'Aures', or 'von Bismarck'
                 'soundFieldType',  'free', ...       % 'free', 'diffuse'
                 'timeVarying',     false);           % false, true
                           
[input,sArgs] = ita_parse_arguments(sArgs,varargin); 

%% Checks

% Matlab version
matlabVer = version('-release');
if str2double(matlabVer(1:4))< 2020
    error('ITA_CALCULATE_SHARPNESS: A newer version of Matlab (at least 2020a) is required for running this function. Use ita_loudness alternatively.')
end

%% Core

data_in = input.timeData;  
fs = input.samplingRate;

if nargout == 0
    if sArgs.timeVarying == true
        if size(data_in,2)>2
            data_calc = data_in(:,1:2);
            disp('=== Displaying only first two channels. ===')
        else
            data_calc = data_in;
        end
        acousticSharpness(data_calc,fs,1,'SoundField',sArgs.soundFieldType,'Weighting',sArgs.weighting,...
                    'TimeVarying', sArgs.timeVarying);%, 'TimeResolution',sArgs.timeResolution);
    else
        disp('No plot possible (per default) for stationary sharpness.')
    end
else
    if size(data_in,2)<=2
        % use input data directly
        [S] = acousticSharpness(data_in,fs,1,'SoundField',sArgs.soundFieldType,'Weighting',sArgs.weighting,...
                    'TimeVarying', sArgs.timeVarying);%, 'TimeResolution',sArgs.timeResolution);
    else
        S = [];
        for idxAudioCh = 1:size(data_in,2)
            [S_temp] = acousticSharpness(data_in(:,idxAudioCh),fs,1,'SoundField',sArgs.soundFieldType,'Weighting',sArgs.weighting,...
                    'TimeVarying', sArgs.timeVarying);%, 'TimeResolution',sArgs.timeResolution);
            S(:,idxAudioCh) = S_temp;
        end
    end
end

%% Set Output

if nargout > 0
%     for idxTime = 1: size(S,1)
%         for idxAudioCh = 1: size(S,2)
%             sharpness(idxTime,idxAudioCh) = itaValue(S(:,idxAudioCh),'sone');
%         end
%     end
    sharpness = S;
    varargout(1) = {sharpness};
end


end