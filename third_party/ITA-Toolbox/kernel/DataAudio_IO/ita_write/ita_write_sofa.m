function result = ita_write_sofa(varargin)
%ITA_WRITE_SOFA - +++ Writes itaObject to Sofa Format +++
%  This function is currently used to save itaHRTF to SOFA format. 
%  Updated to SOFA Version 2.0 -- AES69-2020 (SOFA 2.0)
%  -- www.sofaconventions.org
%  
%  Syntax:
%   ita_write_sofa(hrtfObj,fileName,options)
%
%   options:
%       'dataType' ('') : for later usage to allow spcification of dataryoe
%       when not clear from itaClass
%
%  Example:
%   audioObjOut = ita_write_sofa(hrtfObj,'testHRTF.sofa')
%
%  See also:
%   ita_toolbox_gui, ita_read, ita_write, ita_generate, ita_write_sofa_hrtf
%
%   Reference page in Help browser 
%        <a href="matlab:doc ita_write_sofa">doc ita_write_sofa</a>

% <ITA-Toolbox>
% This file is part of the ITA-Toolbox. Some rights reserved. 
% You can find the license for this m-file in the license.txt file in the ITA-Toolbox folder. 
% </ITA-Toolbox>


% Author: Jan Gerrit Richter -- Email: jan.richter@akustik.rwth-aachen.de
% Created:  30-Sep-2014 
% Update as container for different write function:
% Hark Braren -- hark.braren@akustik.rwth-aachen.de


if nargin == 0 % Return possible argument layout
    result{1}.extension = '*.sofa';
    result{1}.comment = 'SOFA Files (*.sofa)';
    return;
end

sArgs = struct('pos1_data','itaAudio','pos2_filename','char','dataType','','userData','');
[data, filename, sArgs] = ita_parse_arguments(sArgs,varargin); 

%% check if sofa is installed
if ~exist('SOFAstart.m','file')
    error('SOFA not installed. Run ita_sofa_install');
end

if isempty(sArgs.dataType)
    %when not stated explicitly -> derive from data type
    switch class(data)
        case 'itaHRTF'
            sArgs.dataType = 'HRTF';
        otherwise
            error('Unable to determine Sofa type from datatyoe. Check if write function is implemented')
    end
end
    


%%
switch(sArgs.dataType)
    case 'HRTF'
        ita_write_sofa_hrtf(data,filename,'userData',sArgs.userData);
    
%     case 'Directivity'
%         sofaObj = SOFAgetConventions('GeneralTF');   
        
%      case 'SingleRoomDRIR'
%         sofaObj = SOFAgetConventions('SingleRoomDRIR');   

    otherwise
        error('ITA_WRITE_SOFA: Only HRTF Type is defined');
end

result = 1;
end
