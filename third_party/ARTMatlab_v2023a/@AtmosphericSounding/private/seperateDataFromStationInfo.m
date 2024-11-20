function [dataTxt, stationInfoTxt] = seperateDataFromStationInfo(txt)
%The atmospheric sounding data files contain of two parts. First a table
%with the weather data and second the data for the station where these
%measurments where done.
%This function separates both parts and returns them.
%
%Input:
%txt:       Atmospheric sounding data (Cellarray with one string per line)
%
%Outputs:
%dataTxt:           Weather data in text form (Cellarray with one string per line)
%stationInfoTxt:    Station data in text form (Cellarray with one string per line)

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

idxEmptyLines = findEmptyLines(txt);
idxSeparator = find( contains(txt, '--------') );
assert(numel(idxEmptyLines)>2, 'Unsupported file format. Expected at least 3 blank lines.')
assert(numel(idxSeparator)==2, 'Unsupported file format. Expected exactly 2 lines with separators.')

idxDataStart = idxSeparator(2)+1;
idxFirstEmptyLineAfterDataStart = find(idxEmptyLines>idxDataStart, 1, 'first');
idxDataEnd = idxEmptyLines(idxFirstEmptyLineAfterDataStart)-1;

idxStationInfoStart = find( contains(txt, 'Station number') );
idxStationInfoEnd = find( contains(txt, 'Precipitable water') );
if isempty(idxStationInfoEnd)
    idxStationInfoEnd = length(txt);
end

dataTxt = txt(idxDataStart:idxDataEnd);
dataTxt = removeIncompleteLine(dataTxt);

stationInfoTxt = txt(idxStationInfoStart:idxStationInfoEnd);

function dataTxt = removeIncompleteLine(dataTxt)

idxIncomplete = getIdxIncompleteLines(dataTxt);
dataTxt(idxIncomplete) = [];

function idxIncomplete = getIdxIncompleteLines(dataTxt)

idxIncomplete = false(size(dataTxt));
for idx = 1:numel(dataTxt)
    line = dataTxt{idx};
    jumps = diff(line == ' ');
    jPos = sum(jumps == 1);
    jNeg = sum(jumps == -1);
    idxIncomplete(idx) = jNeg ~= 10 || jPos ~=10;
end