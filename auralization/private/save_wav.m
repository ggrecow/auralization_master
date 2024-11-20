function save_wav(inputSignal, fs, AttenuationdB, fileTag, savePath)
% function save_wav(inputSignal, fs, AttenuationdB, fileTag, savePath)
%
% Write .wav file from a given <inputSignal> with a given sampling
% frequency. 
%
% An attenuation factor <AttenuationdB> can be applied 
% in order to avoid clipping. When this is done, the dBFS of the .wav file do no
% correspond anymore to Pascal values (i.e. +1/-1 are not 94 dB SPL),
% so this value needs to be known in order to calibrate the signal for later analysis
%
% INPUTS
%
%       inputSignal : vector
%       signal (Pa) to be writen in .wav format 
%
%       fs (Hz) : scalar
%       sampling rate (Hz) of <inputSignal>
%       
%       AttenuationdB : scalar
%       Attenuation value, in dB SPL, to be applied to the signal. 
%       Ex: If you want to attenuate the signal by 10 dB SPL, 
%       then <AttenuationdB> = -10. In this case, the +1/-1 amplitude 
%       of the written .wav file corresponds to 94-(AttenuationdB) = 104 dBFS.  
%
%       fileTag : string
%       name of the .wav file
%       If <fileTag> is empty, nothing is saved
%
%       savePath : string
%       path where the .wav file is saved
%       
%   OUTPUTS   
%   no outputs are provided in the matlab environment, only the wav written on the desired path  
%   
% Author: Gil Felix Greco, Braunschweig 18.01.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isempty(savePath) % if <save_path> is empty, dont save anything
else
    
    %% attenuate signal
    
%     AttenuationdB = 0; % -10 -> so +1/-1 amplitude corrsponds to 104 dB SPL
    
    % string with the dBFS value according to <AttenuationdB>
    AttenuationTag = ['_FS' sprintf('%.1d',94-AttenuationdB) 'dBSPL']; 
    
    % attenuation factor
    Attenuation = 10^(AttenuationdB/20);    
    
    %% attenuate and write wav file 
    
    % apply attenuation factor (if <AttenuationdB>=0, then <Attenuation>=1
    % and nothing happens, i.e. dBFS=94 or 1 Pa)    
    outputSignal_Attenuated = inputSignal*Attenuation;
    
    % string with name of the file
    fileName = [fileTag AttenuationTag '.wav'];
    
    % write .wav
    audiowrite( [savePath fileName], outputSignal_Attenuated, fs, 'BitsPerSample', 32);
    
end
