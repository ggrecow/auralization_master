function test_ita_psychoacoustic_loudness

sine = itaAudio;
sine.samplingRate = 44100; % standard sampling frequency
sine.trackLength = 10; %length of five seconds
t_sin = sine.timeVector;
sine.time = sin(2*pi* 1000 * t_sin); % frequency is 1kHz
spl = itaValue(40,'dB');
p_0 = itaValue(2*10^(-5),'Pa');
p = p_0 * 10^(spl/20);
sine = sine*p*sqrt(2);

sine2= ita_merge(sine, sine);

% ita_psychoacoustic_loudness(sine2, 'method', 'Zwicker', 'timeVarying', true, 'timeResolution', 'high');

[loudness,itaResultObj_NS] = ita_psychoacoustic_loudness(sine2, 'method', 'Zwicker', 'timeVarying', true, 'timeResolution', 'high');
if isscalar(loudness) 
    X = ['The loudness is calculated to ', num2str(loudness)];
    disp(X)
else 
    X = ['The time-varying loudness is calculated to ', mat2str(loudness), ' sone.'];
    disp(X)
end 
itaResultObj_NS.plot_specific