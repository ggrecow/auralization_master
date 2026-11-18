function results = validate_doppler_delay_line(show)
% VALIDATE_DOPPLER_DELAY_LINE  Verification of <get_doppler_mapping> and
% <apply_doppler_delay_line> against analytical solutions.
%
%   results = validate_doppler_delay_line(show)
%
% Tests (all without toolboxes):
%   1) Flyovers of a point source in a homogeneous medium at rest (c = 340 m/s),
%      propagation times given every 0.1 s, as from the ray tracing: level flight
%      (85 m/s, M = 0.25) at 300 m and at 60 m (fast change of the Doppler factor
%      close to the receiver), and an accelerating, climbing departure. The received
%      signal of a 1 kHz tone is compared sample by sample with the exact solution
%      D^-2 * x(t_e(t_r)), where t_e(t_r) solves t_r = t_e + R(t_e)/c (Newton).
%   2) Same as 1), with random errors added to the propagation times (std 0.1 ms,
%      i.e. 3.4 cm path length) to check the smoothing in <get_doppler_mapping>.
%   3) Doppler factor and convective amplification: D from the propagation times
%      vs. 1 - M*cos(theta), and the level of the received tone vs. -40*log10(D).
%   4) Fractional delay accuracy: constant propagation time with a fractional part of
%      0.5 samples (worst case), level error of tones from 100 Hz to 0.35*fs.
%   5) Stationary source with an integer delay: output must equal the delayed input.
%
% INPUT
%   show : logical, plots results (DEFAULT false)
%
% OUTPUT
%   results : struct with the error measures of each test, and results.pass
%
% Author: Gil Felix Greco, Braunschweig 16.11.2026
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to
%   develop this code under the author's direction. All codes were verified
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1, show = false; end

fs = 48000;
c  = 340;
v  = 85;
M  = v/c;
dt = 0.1;
T  = 20;                         % [s] duration of the emission signal
t0 = 10;                         % [s] closest approach (emission time)
f0 = 1000;                       % [Hz] test tone

fprintf('\n*--------------------------------------------------------------------------*\n');
fprintf('validate_doppler_delay_line: fs = %d Hz, v = %g m/s, c = %g m/s (M = %.3f), dt = %g s\n', fs, v, c, M, dt);

t_e = (0:dt:T).';
n   = (0:round(T*fs)).';
x   = cos( 2*pi*f0*n/fs );       % source signal in emission time

pass = true;

%% tests 1-3: flyovers

% cases: name, altitude h0 [m], speed v0 [m/s], acceleration [m/s^2], climb rate [m/s], tau noise [s]
cases = { 'level flight, h = 300 m',            300, v, 0, 0, 0
          'level flight, h = 300 m, tau noise', 300, v, 0, 0, 1e-4
          'level flight, h = 60 m',             60,  v, 0, 0, 0
          'departure (2 m/s^2, climb 8 m/s)',   200, 80, 2, 8, 0
          'departure, tau noise',               200, 80, 2, 8, 1e-4 };

for k = 1:size(cases,1)
    [name, h0, v0, acc, vz, sig] = cases{k,:};
    traj = @(t) il_trajectory( t, t0, h0, v0, acc, vz );

    [Rk, dRk] = traj( t_e );
    tau  = Rk/c;
    rng(1);
    tauN = tau + sig*randn(size(tau));

    map = get_doppler_mapping( t_e, tauN );

    % output window: inner part (avoids the edges of the smoothing window, 1 s)
    t_out = ( map.t_r(11) : 1/fs : map.t_r(end-10) ).';
    [y, g] = apply_doppler_delay_line( x, fs, map, t_out, 1 );

    % exact solution: t_r = t_e + R(t_e)/c solved with Newton's method
    te_exact = t_out - traj( t_out )/c;
    for it = 1:8
        [R, dR] = traj( te_exact );
        te_exact = te_exact - ( te_exact + R/c - t_out ) ./ ( 1 + dR/c );
    end
    [~, dR] = traj( te_exact );
    D_ex = 1 + dR/c;                            % Doppler factor, = 1 - M*cos(theta)
    y_ex = D_ex.^(-2) .* cos( 2*pi*f0*te_exact );

    err = y - y_ex;
    snr = 10*log10( sum(y_ex.^2) / sum(err.^2) );    % signal-to-error ratio [dB]

    % Doppler factor at the time steps
    Dk = 1 + dRk/c;
    inner = 11:numel(t_e)-10;
    errD = max( abs( map.D(inner) - Dk(inner) ) );

    % level of the tone (convective amplification) in blocks of 50 ms
    L  = round(0.05*fs);
    nb = floor( numel(y)/L );
    Lerr = zeros(nb,1);
    for b = 1:nb
        idx = (b-1)*L + (1:L);
        Lerr(b) = 10*log10( mean(y(idx).^2) / mean(y_ex(idx).^2) );
    end

    results.flyover(k).name = name;
    results.flyover(k).signalToErrorRatio_dB = snr;
    results.flyover(k).maxErrorDopplerFactor = errD;
    results.flyover(k).maxLevelError_dB = max(abs(Lerr));
    results.flyover(k).rangeDopplerFactor = [min(Dk) max(Dk)];
    results.flyover(k).rangeGain_dB = 20*log10([min(g) max(g)]);
    results.flyover(k).tauSmoothingRms = map.tau_smoothing_rms;

    if sig == 0   % waveform must match (timing and phase), not only pitch and level
        ok = snr > 40 && errD < 1e-3 && max(abs(Lerr)) < 0.05;
    else          % random timing errors: pitch and level must be stable
        ok = errD < 5e-3 && max(abs(Lerr)) < 0.1;
    end
    pass = pass && ok;

    fprintf(['- %-36s D = %.3f..%.3f, gain %+.2f..%+.2f dB | signal/error %5.1f dB, ' ...
             'max |D error| %.1e, max level error %.3f dB, smoothing rms %.1e s -> %s\n'], ...
        [name ':'], min(Dk), max(Dk), results.flyover(k).rangeGain_dB, snr, errD, ...
        max(abs(Lerr)), map.tau_smoothing_rms, il_ok(ok));

    if show
        figure('name', ['validate_doppler_delay_line - ' name]);
        subplot(2,1,1);
        plot( t_e, Dk, 'k', t_e, map.D, 'r--' ); grid on;
        xlabel('Emission time, t_e (s)'); ylabel('Doppler factor D (-)');
        legend('exact: 1 - M cos\theta', 'from propagation times', 'Location', 'SE');
        title(name);
        subplot(2,1,2);
        plot( t_out, 20*log10(abs(err)+eps) ); grid on;
        xlabel('Reception time, t_r (s)'); ylabel('|error| (dB re 1)');
    end
end

%% test 4: fractional delay accuracy (constant propagation time, 0.5 samples)

tau0 = 1 + 0.5/fs;                      % [s]
mapC = get_doppler_mapping( t_e, tau0*ones(size(t_e)) );
t_out = ( 2 : 1/fs : T - 1 ).';
fTest = [100 1000 5000 10000 12500 15000 0.35*fs];
Lerr4 = zeros(size(fTest));
for i = 1:numel(fTest)
    xi = cos( 2*pi*fTest(i)*n/fs );
    yi = apply_doppler_delay_line( xi, fs, mapC, t_out, 0 );
    yref = cos( 2*pi*fTest(i)*(t_out - tau0) );
    Lerr4(i) = 10*log10( mean(yi.^2) / mean(yref.^2) );
end
ok4 = max(abs(Lerr4)) < 0.05;
pass = pass && ok4;
results.fractionalDelay.frequencies = fTest;
results.fractionalDelay.levelError_dB = Lerr4;
fprintf('- fractional delay 0.5 samples, tones %g..%g Hz: max level error %.4f dB -> %s\n', ...
    fTest(1), fTest(end), max(abs(Lerr4)), il_ok(ok4));

%% test 5: integer delay (stationary source)

nDelay = 4800;
mapI = get_doppler_mapping( t_e, nDelay/fs*ones(size(t_e)) );
rng(2);
xr = randn(size(n));
t_out = ( 2 : 1/fs : T - 1 ).';
yI = apply_doppler_delay_line( xr, fs, mapI, t_out, 1 );
idx = round( (t_out - nDelay/fs)*fs ) + 1;
errI = max( abs( yI - xr(idx) ) );
ok5 = errI < 1e-6;
pass = pass && ok5;
results.integerDelay.maxAbsError = errI;
fprintf('- integer delay (stationary source): max |error| %.1e -> %s\n', errI, il_ok(ok5));

results.pass = pass;
fprintf('validate_doppler_delay_line: %s\n', il_ok(pass));
fprintf('*--------------------------------------------------------------------------*\n');

end

function [R, dR] = il_trajectory(t, t0, h0, v0, acc, vz)
% distance source-receiver R and its time derivative dR/dt. Receiver at the origin,
% flight along +x, closest approach at t = t0 (x = 0, z = h0)
u  = t - t0;
xs = v0*u + 0.5*acc*u.^2;   vx = v0 + acc*u;
zs = h0 + vz*u;             vzs = vz*ones(size(u));
R  = sqrt( xs.^2 + zs.^2 );
dR = ( xs.*vx + zs.*vzs ) ./ R;
end

function s = il_ok(b)
if b, s = 'OK'; else, s = 'FAILED'; end
end
