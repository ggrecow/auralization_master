function [outData] = lowPassFilter(inData,samplingRate,filterType,passBandFreq,stopBandFreq,maxRipple,stopbandAttenuation)
%LOWPASSFILTER low pass filters inData

[nSamples,nDim] = size(inData);

if nargin < 3
    filterType = 'lowpassfir';
    passBandFreq = 1;
    stopBandFreq = 3.5;
    maxRipple = 1e-4;    
end

if nargin < 7
    stopbandAttenuation = 60;
end
%% extend signal to avoid edge effects
nExtend = max(floor(nSamples/10),200);
extendStart = inData(1,:).*ones(nExtend,nDim);
extendEnd = inData(end,:).*ones(nExtend,nDim);

inData_ext = [extendStart;inData;extendEnd];
%% filter signal

lpFilter = designfilt(filterType,'PassbandFrequency',passBandFreq,'StopbandFrequency',stopBandFreq,...
    'PassbandRipple',maxRipple,'StopbandAttenuation',stopbandAttenuation,'SampleRate',samplingRate);


outData = filter(lpFilter,inData_ext);
outData = outData(nExtend+1:nExtend+nSamples,:);

end

