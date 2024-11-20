function idxEmptyLines = findEmptyLines(textIn)
%Identifies and returns the indices of empty lines of the given string cell.
%
%Input:
%textIn:        Cellarray with one string per text line
%
%Output:
%idxEmptyLines: Indices of empty lines

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

idxEmptyLines = [];
for k = 1:length(textIn)
   line = textIn{k};
   if isempty(line)
       idxEmptyLines = [idxEmptyLines k];
   end
end