function varargout = ita_psychoacoustic_loudness(varargin)
%ITA_CALCULATE_LOUDNESS - Calculates mono, stereo or binaural loudness
%  This is a wrapper function for acousticLoudness.m, which determines 
%  total loudness (N) and specific loudness (N') values for one or two
%  audio channels. It applies 'ISO 532-1' after Zwicker (DIN45631) for mono or stereo 
%  calculations (loudness per channel), or the 'Moore-Glasberg' approach
%  ('ISO 532-2') for binaural loudness. More than two channels are admitted, 
%  monaural calculations are then performed. 
%  Requirement: Matlab 2020a or later.
%
%  Syntax:
%   [double_N,itaResultObj_NS] = ita_calculate_loudness(audioObjIn, options)
%
%   Options (default):
%           'method' ('Zwicker'/'ISO 532-1') : loudness calculation method ('Moore-Glasberg'/'ISO 532-2')
%           'soundFieldType' ('free')        : field conditions ('diffuse', 'eardrum', 'earphones')
%           'timeVarying' (false)            : time-dependency (true)
%           'timeResolution' ('standard')    : time window size ('standard': 2 ms, 'high': 0.5 ms) 
% 
%  Example:
%   [N_tot, N_spec] = ita_calculate_loudness(audioObjIn, 'soundFieldType','free', 'method','Zwicker','timeVarying',true)
%
%  See also:
%   ita_toolbox_gui, ita_read, ita_write, ita_generate
%
%   Reference page in Help browser 
%        <a href="matlab:doc ita_loudness2">doc ita_loudness2</a>

% <ITA-Toolbox>
% This file is part of the ITA-Toolbox. Some rights reserved. 
% You can find the license for this m-file in the license.txt file in the ITA-Toolbox folder. 
% </ITA-Toolbox>


% Author: ShaimaaDoma -- Email: sdo@akustik.rwth-aachen.de
% Created:  04-Oct-2021 

% TODO: option to hand over 1/3 octave BP-filtered signal (28 or 29 bands)
% -- needed?

%% Initialization and Input Parsing
sArgs   = struct('pos1_data',       'itaAudio', ...
                 'method',          'Zwicker', ...    % 'Zwicker' ('ISO 532-1'), 'Moore-Glasberg' ('ISO 532-2') 
                 'soundFieldType',  'free', ...       % 'free', 'diffuse', 'eardrum'* or 'earphones'*
                 'timeVarying',     false,...         % false, true**
                 'timeResolution',  'standard');      % 'standard': 2 ms, 'high':0.5 ms
[input,sArgs] = ita_parse_arguments(sArgs,varargin); 

% Note:  * only feasible for Moore-Glasberg method
%       ** only feasible for Zwicker method


%% Checks
% Matlab version
matlabVer = version('-release');
if str2double(matlabVer(1:4))< 2020
    error('ITA_CALCULATE_LOUDNESS: A newer version of Matlab (at least 2020a) is required for running this function. You can use ita_loudness alternatively.')
end

% Unit
if ~all(strcmp(input.channelUnits(:), 'Pa'))
    ita_verbose_info('Input itaAudio object is assumed to have Pa as unit. Double check please, a scaling factor might be missing.',0)
    %warning([thisFuncStr,]);
end
%% Core

switch sArgs.method
    case {'Zwicker','ISO 532-1'} % monaural model (for use on mono or stereo data)
        method = 'ISO 532-1'; 
        barkVector = (0.1:0.1:24)';
    case {'Moore-Glasberg','ISO 532-2'} % binaural model
        method = 'ISO 532-2';
        barkVector = (1.8:0.1:38.9)';
    otherwise
        error([thisFuncStr,'Unknown loudness calculation method. Please choose between ''Zwicker'' and  ''Moore-Glasberg''.']);      
end

fs = input.samplingRate;
data_in = input.timeData;  

if nargout == 0
    if size(data_in,2)>2
        for idxAudioCh = 1:size(data_in,2)
            figure
           acousticLoudness(data_in(:,idxAudioCh),fs,1,'Method',method,'SoundField',sArgs.soundFieldType,...
                    'TimeVarying', sArgs.timeVarying, 'TimeResolution',sArgs.timeResolution);
        end
    else
        acousticLoudness(data_in,fs,1,'Method',method,'SoundField',sArgs.soundFieldType,...
                    'TimeVarying', sArgs.timeVarying, 'TimeResolution',sArgs.timeResolution);
    end
 
else
    if size(data_in,2)<=2
        % use input data directly -- interpretation as monaural or binaural possible
        [N, NS] = acousticLoudness(data_in,fs,1,'Method',method,'SoundField',sArgs.soundFieldType,...
                    'TimeVarying', sArgs.timeVarying, 'TimeResolution',sArgs.timeResolution);
    else
        % assume multiple monaural calculations
        N = []; NS = [];
        for idxAudioCh = 1:size(data_in,2)
            [N_temp, NS_temp] = acousticLoudness(data_in(:,idxAudioCh),fs,1,'Method',method,'SoundField',sArgs.soundFieldType,...
                    'TimeVarying', sArgs.timeVarying, 'TimeResolution',sArgs.timeResolution);
            N(:,idxAudioCh) = N_temp;
            NS(:,:,idxAudioCh) = NS_temp;
        end
    end
end

%% Set Output
if nargout > 0
    totalLoudness = N;
    varargout(1) = {totalLoudness}; 
if nargout > 1
    % freqVector         = (0.1: 0.1: 24)';
    if strcmp(sArgs.timeResolution, 'standard')
        time_res = 2e-3;
    elseif strcmp(sArgs.timeResolution, 'high')
        time_res = 5e-4;
    else
        ita_verbose_info('Something might have changed in the time resolution parameters of the Matlab function acousticLoudness, please check!')
    end

    specLoudness = itaSpecificResult(NS, 0:time_res:(numel(NS(:,1))-1)*time_res, barkVector);
    for idxAudioCh = 1:size(data_in,2)
        specLoudness.channelUnits(idxAudioCh,:)    = {'sone/Bark'};
        specLoudness.channelDataType(idxAudioCh,:)    = {'Specific Loudness'};
    end
    varargout(2) = {specLoudness};
end
end


end

% switch sArgs.method
%     case {'Zwicker','ISO 532-1'} % monaural model
%         method = 'ISO 532-1';
%         x = input.timeData;
%     case {'Moore-Glasberg','ISO 532-2'} % binaural model
%         method = 'ISO 532-2';
%         if isa(input, 'itaHRTF')
%             % Make sure that left and right ear data are arranged in an alternating order
%             x1 = input.getEar('L').timeData;
%             x2 = input.getEar('R').timeData;
%             x = [x1;x2]; x = reshape(x,size(x1,1),2*size(x1,2));
%         else
%             x = input.timeData; 
%         end
%     otherwise
%         error([thisFuncStr,'Unknown loudness calculation method. Please choose between ''Zwicker'' and  ''Moore-Glasberg''.']);      
% end
