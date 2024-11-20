function [f1, fm, f2] = get_octave_bands(b)
% function [f1, fm, f2] = get_octave_bands(b)
%
% Get 1/b octave bands according to:
% 
% [1] DIN 61260-1:2014
%
%   INPUT
%
%   b : scalar
%   Bandwidth metric, e.g., b=1 for octave filter or b=3 third octave filter.
%
%   OUTPUTS
%   fm - central frequencies of each 1/3 octave band
%   f1 - lower frequencies of each 1/3 octave band
%   f2 - upper frequencies of each 1/3 octave band
%
% Gil Felix Greco, Braunschweig 18.04.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G   = 10^(3/10);    % Eq. 1 in [1]
fr  = 1000; % reference frequency of 1 kHz
% b   = 3;    % Bandwidth metric, e.g., b=1 for octave filter or b=3 third octave filter.

% x = -16:13; % from Tabelle E.1 - gives bands from 25 Hz to 20 kHz
x = -18:14; % gives bands from 16 Hz to 25 kHz

% exact central freqs.
fm = fr.*G.^(x/b); % Eq. 2 in [1] 

% lower freqs.
f1 = fm.*G.^-(1 ./(2.*b));  % Eq. 4 in [1] 

% upper freqs.
f2 = fm.*G.^(1 ./(2.*b));   % Eq. 5 in [1] 

% fnom = [ 25 31.5 40, 50 63 80, 100 125 160, 200 250 315, 400 500 630, ...
%     800 1000 1250, 1600 2000 2500, 3150 4000 5000, 6300 8000 10000, ...
%     12500 16000 20000 ];    % nominal center freq - preferred for freq. labeling


end