function [outTF_sum,outTF_separated] = applyAtmoAttenuation(urbanTF,attenuation,freqVector)
%APPLYATMOATTENUATION applies atmospheric TF parameters to urban TF


if urbanTF.nBins == length(freqVector)
    %% apply attenuation directly
    outTF_separated = itaAudio();
    
    outTF_separated.freqData = urbanTF.freqData .* attenuation';
    outTF_sum = outTF_separated.sum;
    
else    
    %% apply third octave band filter and multiply complex attenuation factor

    freqRange = [freqVector(1), freqVector(length(freqVector))];

    for idChannel = 1:urbanTF.nChannels

        filteredTF = ita_fractional_octavebands(urbanTF.ch(idChannel),'bandsperoctave',3,'freqRange',freqRange);
        filteredTF = filteredTF .* attenuation(idChannel,:);

        if idChannel == 1
            outTF_sum = filteredTF.sum;
            outTF_separated = filteredTF.sum;
        else
            outTF_sum = ita_add(outTF_sum,filteredTF.sum);
            outTF_separated = ita_merge(outTF_separated,filteredTF.sum);
        end

        outTF_separated.channelNames{idChannel} = ['Path ',num2str(idChannel)];

    end
    % filteredTF = ita_fractional_octavebands(urbanTF,'bandsperoctave',3,'freqRange',freqRange);
    outTF_sum.channelNames{1} = 'Combined Paths';

end
end

