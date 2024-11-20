function test_ita_sharpness()

% <ITA-Toolbox>
% This file is part of the application Psychoacoustics for the ITA-Toolbox. All rights reserved. 
% You can find the license for this m-file in the application folder. 
% </ITA-Toolbox>


%% generate standard test signal
flatCBlevel = ita_generate('pinknoise',18.5,44100,18);
flatCBlevel.channelUnits = {'Pa'};

%% edit standard signal
schmalbandRauschen = flatCBlevel;
schmalbandRauschen.freqData([1:schmalbandRauschen.freq2index(920) schmalbandRauschen.freq2index(1080):end]) = 0;


% terzen = ita_spk2frequencybands(ita_fft(schmalbandRauschen), 'fraction',3,'method','added','limits',[25 12500]);
% terzen.bar


%% evaluate sharpness
S = ita_sharpness(schmalbandRauschen);
 X = ['The sharpness is calculated to ', num2str(S)];
disp(X), 

%% evaluate sharpness calculated with test_ita_calculate_sharpness
S = ita_calculate_sharpness(schmalbandRauschen);
Y = ['The sharpness is calculated to ', num2str(S), ' acum.'];
disp(Y), 

% ita_sharpness(flatCBlevel)