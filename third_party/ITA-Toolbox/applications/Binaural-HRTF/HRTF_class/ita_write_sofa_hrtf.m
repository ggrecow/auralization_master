function varargout = ita_write_sofa_hrtf(varargin)
%ITA_WRITE_SOFA_HRTF - +++ Write HRTF as SOFA file +++
%  This function writes a itaHRTF into a .sofa file following the SOFA 
%
%  Syntax:
%   ita_write_sofa_hrtf(hrtfObj, fileName)
%
%  Options:
%    userData    ([]): struct with user specific SOFA entries see
%                      <userDataFields> below      
%                      e.g. userData = struct('GLOBAL_AuthorContact','Max Mustermann: mmu@xyz.de',...
%                                             'GLOBAL_Comment','Example Comment')


%  See also:
%   ita_read_sofa_hrtf, ita_write, itaHRTF
%
%   Reference page in Help browser 
%        <a href="matlab:doc ita_write_sofa_hrtf">doc ita_write_sofa_hrtf</a>

% <ITA-Toolbox>
% This file is part of the ITA-Toolbox. Some rights reserved. 
% You can find the license for this m-file in the license.txt file in the ITA-Toolbox folder. 
% </ITA-Toolbox>


% Author: Hark Braren  -- Email: hark.braren@akustik.rwth-aachen.de
% Created:  12-Apr-2021 

userDataFields = {'GLOBAL_AuthorContact',... %Who created it, automatically takes info from toolbox preferences
                  'GLOBAL_Comment',...
                  'GLOBAL_History',...
                  'GLOBAL_License',...       % e.g. CC 4.0 BY - http://creativecommons.org/licenses/by/4.0/
                  'GLOBAL_ListenerShortName',...  % e.g. KEMAR
                  'GLOBAL_Organization',...  % (IHTA)
                  'GLOBAL_References',...    % corresonding paper
                  'GLOBAL_DateCreated',...
                  'GLOBAL_DateModified',...
                  'GLOBAL_Title',...
                  'GLOBAL_DatabaseName',...  %
                  'GLOBAL_Origin',...
                  'GLOBAL_RoomType',...      %(freefield)
                  };

%% Initialization and Input Parsing
if nargin == 0 % Return possible argument layout
    result{1}.extension = '*.sofa';
    result{1}.comment = 'SOFA Files (*.sofa)';
    return;
end

sArgs = struct('pos1_data','itaHRTF','pos2_filename','char','userData',[],'delaySamples',[0,0]);
[data, fileName, sArgs] = ita_parse_arguments(sArgs,varargin); 

%% check if sofa is installed
if ~exist('SOFAstart.m','file')
    error('SOFA not installed. Run ita_sofa_install');
end

SOFAstart

%% get sofa definitions and domain specific data
switch data.domain
    case 'time'
        Obj=SOFAgetConventions('SimpleFreeFieldHRIR');
        
        %time data
        Obj.Data.IR = zeros(data.nDirections, 2, data.nSamples);
        Obj.Data.IR(:,1,:) = data.getEar('L').timeData.';
        Obj.Data.IR(:,2,:) = data.getEar('R').timeData.';            
        
        %data info
        Obj.Data.SamplingRate = data.samplingRate;
        Obj.Data.SamplingRate_Units = 'hertz';
        
    case 'freq'
        Obj=SOFAgetConventions('SimpleFreeFieldHRTF');
        
        %freq data 
        Obj.Data.Real = zeros(data.nDirections, 2, data.nBins);
        Obj.Data.Imag = zeros(data.nDirections, 2, data.nBins);
        
        Obj.Data.Real(:,1,:) = real(data.getEar('L').freqData).';
        Obj.Data.Real(:,2,:) = real(data.getEar('R').freqData).';
        Obj.Data.Imag(:,1,:) = imag(data.getEar('L').freqData).';
        Obj.Data.Imag(:,2,:) = imag(data.getEar('R').freqData).';
        
        %data info
        Obj.N = data.freqVector;
        Obj.N_LongName = 'frequency';
        Obj.N_Units = 'hertz';
end
   

%% General Info
Obj.GLOBAL_ApplicationName = 'ITA-Toolbox';
Obj.GLOBAL_Organization    = 'Institute for Hearing Technology and Acoustics, RWTH Aachen University';
Obj.GLOBAL_ApplicationVersion = num2str(ita_toolbox_version_number);

%The following will be overwritten with userdata is supplied
Obj.GLOBAL_DateCreated   = date; 
Obj.GLOBAL_AuthorContact = [ita_preferences('AuthorStr') ': ' ita_preferences('EmailStr')];


%% HRTF Info
%head position 
Obj.ListenerPosition_Type  = 'cartesian';
Obj.ListenerPosition_Units = 'metre';
Obj.ListenerPosition       = [0 0 0];

% ear position for each measurement [0 +w 0; 0 -w 0] in hrtf.objectCoordinates
Obj.ReceiverPosition_Type  = 'cartesian';
Obj.ReceiverPosition_Units = 'metre';
if ~isempty(data.objectCoordinates.cart)
    Obj.ReceiverPosition       = repmat(data.objectCoordinates.cart,1,1,data.nDirections);
else
    Obj.ReceiverPosition       = repmat([0 0.09 0; 0 -0.09 0],1,1,data.nDirections); %default per SOFA def
end

%head orientation
Obj.ListenerView_Type  = 'cartesian';
Obj.ListenerView_Units = 'metre';
if ~isempty(data.objectViewVector.cart)
    Obj.ListenerView       = data.objectViewVector.cart;
end
if ~isempty(data.objectUpVector.cart)
    Obj.ListenerUp         = data.objectUpVector.cart;
end

%Loudspeaker positions in head-related coordinates
Obj.SourcePosition_Type  = 'spherical';
Obj.SourcePosition_Units = 'degree, degree, metre';
Obj.SourcePosition       = sofaSphericalCoordiantes(data.dirCoord);

if all(size(sArgs.delaySamples) == [1, 1])
    Obj.Data.Delay = repmat(sArgs.delaySamples,1,2);
elseif all(size(sArgs.delaySamples) == [1, 2])
    Obj.Data.Delay = sArgs.delaySamples;
else
    error('Unknown size of Delay Samples, should be single integer, or 1x2 array');
end
%% add Userdata
if ~isempty(sArgs.userData)
    for iField = 1:numel(userDataFields)
        if ismember(userDataFields{iField},fieldnames(sArgs.userData))
            Obj.(userDataFields{iField}) = sArgs.userData.(userDataFields{iField});
        end
    end
end

%% check common errors
if ~strcmp(Obj.GLOBAL_RoomType,'free field')
    ita_verbose_info('For HRTFs, ''GLOBAL_RoomType'' must be set to ''free field''. I''ll do that for you now',1)
    Obj.GLOBAL_RoomType = 'free field';
end
    

%% update and write obj
Obj=SOFAupdateDimensions(Obj);
SOFAsave(fileName,Obj);
end

function sofaCoord = sofaSphericalCoordiantes(coords)
%transform from zenith to elevation angle

    r = coords.r;
    elevation = 90-coords.theta_deg;
    azimuth = coords.phi_deg;
    
    sofaCoord = [azimuth,elevation,r];
end
