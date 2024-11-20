function text = readTextLinewise(filename)
%Reads a textfile linewise.
%
%Input:
%filename:  Name of text file
%
%Output:
%text:      Cellarray with one string per text line

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

fileID = fopen(filename);
text = textscan(fileID,'%s','Delimiter','\n');
text = text{1};
fclose(fileID);