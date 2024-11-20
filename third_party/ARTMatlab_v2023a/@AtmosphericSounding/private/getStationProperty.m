function strStationProp = getStationProperty(stationInfoTxt, propID)
%Returns a station property giving the corresponding data in text form.

% -------------------------------------------------------------------------
%                        ____  __________  _______
%                       //  / //__   ___/ //  _   |
%                      //  /    //  /    //  /_|  |
%                     //  /    //  /    //  ___   |
%                    //__/    //__/    //__/   |__|
%
% -------------------------------------------------------------------------
%                  ARTMatlab - ITAGeometricalAcoustics
%        (c) Copyright Institute of Technical Acoustics (ITA)
%  This file is part of the ARTMatlab application. Some rights reserved.
% You can find the license in the LICENSE.md file in the ARTMatlab folder.
%--------------------------------------------------------------------------

strStationProp = '';
for idxLine = 1:length(stationInfoTxt)
    line = stationInfoTxt{idxLine};
    splitLine = strsplit(line, ': ');
    propId = splitLine{1};
    propVal = splitLine{2};
    
    if(strcmp(propId,propID))
        strStationProp = propVal;
        break;
    end
end