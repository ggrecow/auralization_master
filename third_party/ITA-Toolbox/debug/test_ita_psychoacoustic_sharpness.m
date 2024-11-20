function test_ita_psychoacoustic_sharpness

%% create reference signal
% narrow band random noise bandfiltered for 920 - 1080 Hz
noise = ita_generate('noise',1,44100,15);
spl = itaValue(60,'dB');
p_0 = itaValue(2*10^(-5),'Pa');
p = p_0 * 10^(spl/20);
noise = noise*p*sqrt(2);
f_low  = 920;
f_high = 1080;
noise_filtered = ita_filter_bandpass(noise,'lower',f_low,'upper',f_high);
noise_filtered.channelNames{1} = 'Sharpness Reference Signal';

% S = ita_psychoacoustic_sharpness(noise_filtered, 'weighting', 'von Bismarck');
S = ita_psychoacoustic_sharpness(noise_filtered,'soundFieldType', 'free', 'timeVarying', false);

if isscalar(S) 
    X = ['The sharpness is calculated to ', num2str(S), ' acum.'];
disp(X), 
else 
    X = ['The time-varying sharpness is calculated to ', mat2str(S), ' acum.'];
    disp(X)

end 


end 

