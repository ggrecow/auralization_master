function [idxLines, placeInLine] = searchInText(textIn, search)
%Searches a string within a cellarray of one string per textline and
%returns the indices of lines and for each line the index within this line.
%
%Inputs:
%textIn:    Cellarray of one string per textline
%search:    String pattern which is searched
%
%Outputs:
%idxLines:      Indices of lines which contain the given pattern
%placeInLine:   Indices of the place of the pattern within each line

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

idxLines = [];
placeInLine = [];
for k = 1:length(textIn)
    line = textIn{k};
    idxCurrentLine = strfind(line, search);
    if ~isempty(idxCurrentLine)
        idxLines = [idxLines, k];
        placeInLine = [placeInLine idxCurrentLine(1)];
    end
end