function test_ita_psychoacoustic_fluctuationStrength

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
cosine.time = cos(2*pi*4*t_cos);

% total signal
modulated_signal = (1+cosine).*sine* p*sqrt(2);
modulated_signal.channelNames{1} = 'Fluctuation Strength Reference Signal';

modulated_signal2= ita_merge(modulated_signal, modulated_signal, modulated_signal);

[fluct_tot, fluct_spec, f_mod] = ita_psychoacoustic_fluctuationStrength(modulated_signal2, 'modFreq', 4);

X = ['The fluctuation strength is calculated to ', num2str(fluct_tot), ' vacil.'];
disp(X)

fluct_spec.plot_specific


end 

