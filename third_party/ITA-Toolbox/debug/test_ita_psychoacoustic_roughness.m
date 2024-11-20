function test_ita_psychoacoustic_roughness

%% create reference signal
% sine with f=1kHz 
sine = itaAudio;
sine.samplingRate = 44100; % standard sampling frequency
sine.trackLength = 10; %length of five seconds
t_sin = sine.timeVector;
sine.time = sin(2*pi* 1000 * t_sin); % frequency is 1kHz
spl = itaValue(60,'dB');
p_0 = itaValue(2*10^(-5),'Pa');
p = p_0 * 10^(spl/20);

% modulation cosine
cosine = itaAudio;
cosine.samplingRate = 44100;
cosine.trackLength = 10;
t_cos = cosine.timeVector;
cosine.time = cos(2*pi*70*t_cos); %modulation frequency of 70Hz

% total signal
modulated_signal = (1+cosine).*sine* p*sqrt(2);
modulated_signal.channelNames{1} = 'Roughness Reference Signal';

modulated_signal2 = ita_merge(modulated_signal, modulated_signal, modulated_signal);

[R_tot, R_spec, f_mod] = ita_psychoacoustic_roughness(modulated_signal2, 'soundFieldType','free');

X = ['The roughness is calculated to ', num2str(R_tot), ' asper.'];
disp(X)

R_spec.plot_specific

end 

