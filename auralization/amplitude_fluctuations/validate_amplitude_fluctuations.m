% VALIDATE_AMPLITUDE_FLUCTUATIONS
% Checks the synthesis of turbulence-induced amplitude fluctuations
% (apply_amplitude_fluctuations, get_amplitude_modulation), in the spirit of
% validate_decorrelation_filter.m:
%
%   1) exact reconstruction for sigma = 0, reproducibility with the same seed
%   2) per 1/3-octave band: standard deviation of the level fluctuations of
%      a long white-noise signal vs. the target sigma(f), and the change of
%      the mean energy (must be ~0 dB)
%   3) modulation spectrum, analyzed as in Lincke & Pieren (2024), Sec. 3.1.2
%      (446-1120 Hz, 100 ms Hann, 15 s Hann detrending, 20 s windows, octave
%      bands 0.25-8 Hz, normalized by the total std), compared with their
%      measured Fig. 11 and Lincke et al. (2025) Fig. 4
%   4) correlation of the level fluctuations across frequency (500 Hz vs 2 kHz)
%   5) interaction with the coherence-loss filters: level fluctuations caused
%      by the decorrelation filters alone (redrawn every 0.1 s), by the
%      amplitude fluctuations alone, and by both, for a direct + reflected
%      signal; and the combined change of the mean energy
%
% The target sigma(f) is that of case III of Lincke et al. (2025) (A320 at
% 700 m, meteo condition C) at overhead, from get_log_amplitude_fluctuation.
%
% EXPECTED OUTCOME (18.11.2026, GNU Octave 8.4, seeds as below)
%   1) reconstruction error 2e-15, identical outputs for identical seeds
%   2) synthesized std ~5 % below the target in all bands (the 100 ms Leq
%      averages the 0.5 Hz modulation slightly; target 1.4-3.5 dB), mean energy
%      change -0.02 to -0.11 dB over 300 s. Energy is preserved in EXPECTATION;
%      a single realization deviates, the more the larger sigma (log-normal gain).
%   3) synthesized 0.45 / 0.47 / 0.46 / 0.36 / 0.26 / 0.18 (0.25-8 Hz octaves);
%      LP24 Fig. 11 (read off the figure, approx.): 0.40 / 0.46 / 0.45 / 0.33 /
%      0.21 / 0.14. LP24 attribute the 4 and 8 Hz values largely to background noise.
%   4) correlation 0.999 (single random process for all frequencies)
%   5) decorrelation filters alone (redrawn every 0.1 s): 0.3-1.9 dB level std
%      (hr = 1.6 and 10 m, 250 Hz - 4 kHz), amplitude fluctuations: 2.7-3.2 dB.
%      Both together increase the std only slightly over AF alone (e.g. 2.95 ->
%      3.22 dB at 500 Hz, hr = 10 m; most at 250 Hz: 2.73 -> 3.40 dB). For weak
%      turbulence (sigma ~ 1 dB) the contribution of the decorrelation filters is
%      comparable to the target. Mean energy change -0.04 to -0.15 dB in 60 s.
%
% Author: Gil Felix Greco, Braunschweig (29.10.2026)
% 18.11.2026: expected outcome updated (saturation of the total in get_log_amplitude_fluctuation)
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to 
%   develop this code under the author's direction. All codes were verified 
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

fs   = 48e3;
Tsig = 300;                                   % duration of the test signal [s]
rs   = RandStream('mt19937ar', 'Seed', 1);
x    = randn(rs, round(Tsig*fs), 1);

% target: LKP25 case III at overhead (meteo C)
atm   = struct('ustar', 0.28, 'QH', 109, 'zi', 1040, 'Ts', 295.15);
fCtrl = 10.^((13:43)/10);                     % 1/3-octave centres, 20 Hz - 20 kHz
sigCtrl = get_log_amplitude_fluctuation(700, 10, 690, atm, fCtrl);

%% 1) reconstruction and reproducibility

y0 = apply_amplitude_fluctuations(x(1:10*fs), fs, 0, fCtrl, zeros(size(fCtrl)));
fprintf('1) reconstruction error (sigma = 0): %.1e\n', max(abs(y0 - x(1:10*fs))));
ya = apply_amplitude_fluctuations(x(1:10*fs), fs, 0, fCtrl, sigCtrl, struct('seed', 7));
yb = apply_amplitude_fluctuations(x(1:10*fs), fs, 0, fCtrl, sigCtrl, struct('seed', 7));
yc = apply_amplitude_fluctuations(x(1:10*fs), fs, 0, fCtrl, sigCtrl, struct('seed', 8));
fprintf('   same seed: max diff %.1e | different seed: max diff %.1e\n\n', ...
    max(abs(ya - yb)), max(abs(ya - yc)));

%% 2) std per 1/3-octave band and energy neutrality

y = apply_amplitude_fluctuations(x, fs, 0, fCtrl, sigCtrl, struct('seed', 1));

fB = 10.^((20:39)/10);                        % 100 Hz - 8 kHz
[Lx, ~, Px] = third_octave_level_history(x, fs, fB, 0.1);
[Ly, ~, Py] = third_octave_level_history(y, fs, fB, 0.1);
dLy    = Ly - Lx;                             % same input frames -> noise of the estimate cancels
sigEst = std(dLy, 0, 1);
sigTgt = interp1(log(fCtrl), sigCtrl, log(fB));
dE     = 10*log10(sum(Py, 1) ./ sum(Px, 1));

fprintf('2) band | target std | synthesized std | mean energy change\n');
for b = 1:numel(fB)
    fprintf('   %6.0f Hz | %5.2f dB | %5.2f dB | %+5.2f dB\n', fB(b), sigTgt(b), sigEst(b), dE(b));
end

figure('Position', [100 100 1200 380])
subplot(1,3,1); hold on; box on; grid on
semilogx(fB, sigTgt, 'k-', 'LineWidth', 1.4, 'DisplayName', 'target');
semilogx(fB, sigEst, 'o', 'LineWidth', 1.2, 'DisplayName', 'synthesized (100 ms L_{eq})');
semilogx(fB, dE, 's', 'LineWidth', 1.2, 'DisplayName', 'mean energy change');
set(gca, 'XScale', 'log'); xlim([80 1e4]);
xlabel('1/3-octave band (Hz)'); ylabel('dB'); legend('Location', 'northwest')
title('std of level fluctuations')

%% 3) modulation spectrum (processing of Lincke & Pieren 2024, Sec. 3.1.2)

fsL  = 100;                                   % rate of the level time history [Hz]
hopL = fs/fsL;
Nfft = 2^nextpow2(numel(y));
fF   = (0:Nfft-1).' * fs/Nfft; fF = min(fF, fs - fF);
Y    = fft(y, Nfft);
Y(fF < 446 | fF > 1120) = 0;                  % 446-1120 Hz band pass (brick-wall)
yb_  = real(ifft(Y)); yb_ = yb_(1:numel(y));

E    = sum(reshape(yb_(1:floor(numel(yb_)/hopL)*hopL).^2, hopL, []), 1).';  % 10 ms energies
wH   = 0.5 - 0.5*cos(2*pi*(1:9).'/10);        % 100 ms Hann at 10 ms resolution
Lp   = 10*log10(conv(E, wH/sum(wH), 'same') + realmin);
w15  = 0.5 - 0.5*cos(2*pi*(1:15*fsL-1).'/(15*fsL));
Fp   = Lp - conv(Lp, w15/sum(w15), 'same');   % detrended level fluctuations
Fp   = Fp(15*fsL:end-15*fsL);                 % discard edges of the moving average

fMod = [0.25 0.5 1 2 4 8];
nWin = 20*fsL;
nSeg = floor(numel(Fp)/nWin);
R    = zeros(nSeg, numel(fMod));
fM   = (0:nWin-1).' * fsL/nWin; fM = min(fM, fsL - fM);
for s = 1:nSeg
    seg = Fp((s-1)*nWin + (1:nWin)); seg = seg - mean(seg);
    S   = fft(seg);
    for b = 1:numel(fMod)
        Sb = S; Sb(fM < fMod(b)/sqrt(2) | fM >= fMod(b)*sqrt(2)) = 0;
        R(s, b) = std(real(ifft(Sb))) / std(seg);
    end
end
modSyn  = mean(R, 1);
modLP24 = [0.40 0.46 0.45 0.33 0.21 0.14];    % LP24 Fig. 11, read off the figure (approx.)

fprintf('\n3) modulation band | synthesized | LP24 Fig. 11 (approx.)\n');
for b = 1:numel(fMod)
    fprintf('   %5.2f Hz | %.2f | %.2f\n', fMod(b), modSyn(b), modLP24(b));
end
fprintf('   (%d windows of 20 s; total std of the detrended level: %.2f dB)\n', nSeg, std(Fp));

subplot(1,3,2); hold on; box on; grid on
plot(1:numel(fMod), modLP24, 'k-o', 'LineWidth', 1.4, 'DisplayName', 'LP24 Fig. 11 (approx.)');
he = errorbar(1:numel(fMod), modSyn, std(R, 0, 1), 's-');
set(he, 'LineWidth', 1.2, 'DisplayName', 'synthesized');
set(gca, 'XTick', 1:numel(fMod), 'XTickLabel', {'0.25','0.5','1','2','4','8'});
xlabel('Modulation frequency (Hz)'); ylabel('normalized std'); ylim([0 0.6]); legend
title('modulation spectrum, 446-1120 Hz')

%% 4) correlation across frequency

i500 = find(fB >= 500, 1); i2k = find(fB >= 2000, 1);
k1s  = ones(10, 1)/10;                         % 1 s moving average (reduces estimation noise)
c    = corrcoef(conv(dLy(:, i500), k1s, 'valid'), conv(dLy(:, i2k), k1s, 'valid'));
fprintf('\n4) correlation of level fluctuations, 500 Hz vs 2 kHz (1 s average): %.3f\n', c(1,2));

%% 5) interaction with the coherence-loss filters (direct + reflected path)

Tint = 60;                                     % [s]
xi   = x(1:Tint*fs);
hs = 700; xr = 1500; c0 = 20.05*sqrt(atm.Ts);
L_block = round(0.1*fs);                       % 0.1 s epochs, as in apply_propagation
fBi = [250 500 1000 2000 4000];

fprintf('\n5) interaction with coherence loss (hs = %d m, xr = %d m, |Q| = 1, meteo C)\n', hs, xr);
fprintf('   hr  | band    | target AF | coh. loss only | AF only | both  | energy change (both)\n');

subplot(1,3,3); hold on; box on; grid on
for hr = [1.6 10]
    R1 = sqrt(xr^2 + (hs - hr)^2);
    R2 = sqrt(xr^2 + (hs + hr)^2);
    alpha2 = atan2(hs + hr, xr);
    w  = get_coherence_factor_slanted(hs, hr, R1, alpha2, atm);
    sigAF = get_log_amplitude_fluctuation(hs, hr, R1, atm, fCtrl);

    d    = round((R2 - R1)/c0*fs);             % relative delay of the reflected path [samples]
    N    = 2^11;                               % decorrelation filter length (get_decorrelation_filter default)
    xdir = [zeros(N/2, 1); xi];                % direct path, delayed by the filter centre
    xref = [zeros(d, 1); xi];                  % reflected path
    nOut = numel(xi) + N + d;
    yCoh = zeros(nOut, 1); yCl = zeros(nOut, 1);
    rsF  = RandStream('mrg32k3a', 'Seed', 1);
    nBl  = floor(numel(xref)/L_block);
    for b = 1:nBl
        idx = (b-1)*L_block + (1:L_block);
        rsF.Substream = b;
        [~, bhat] = get_decorrelation_filter(w, 2^14, fs, rsF);
        seg = conv(xref(idx), bhat);
        yCl(idx(1):idx(1)+numel(seg)-1) = yCl(idx(1):idx(1)+numel(seg)-1) + seg;
    end
    yCoh(N/2 + (1:numel(xref))) = xref;        % reflected path without decorrelation
    yDir = zeros(nOut, 1); yDir(1:numel(xdir)) = xdir;
    yA   = yDir + yCoh;                        % coherent sum (reference)
    yB   = yDir + yCl;                         % coherence loss only
    yC   = apply_amplitude_fluctuations(yA, fs, 0, fCtrl, sigAF, struct('seed', 1));   % AF only
    yD   = apply_amplitude_fluctuations(yB, fs, 0, fCtrl, sigAF, struct('seed', 1));   % both

    LA = third_octave_level_history(yA, fs, fBi, 0.1);
    [LB, ~, ~]  = third_octave_level_history(yB, fs, fBi, 0.1);
    [LC, ~, ~]  = third_octave_level_history(yC, fs, fBi, 0.1);
    [LD, ~, PD] = third_octave_level_history(yD, fs, fBi, 0.1);
    [~, ~, PA]  = third_octave_level_history(yA, fs, fBi, 0.1);
    [~, ~, PB]  = third_octave_level_history(yB, fs, fBi, 0.1);
    sB = std(LB - LA, 0, 1); sC = std(LC - LA, 0, 1); sD = std(LD - LA, 0, 1);
    sT = interp1(log(fCtrl), sigAF, log(fBi));
    eD = 10*log10(sum(PD, 1)./sum(PB, 1));     % AF on top of coherence loss
    for b = 1:numel(fBi)
        fprintf('   %4.1f | %5.0f Hz | %5.2f dB | %5.2f dB | %5.2f dB | %5.2f dB | %+5.2f dB\n', ...
            hr, fBi(b), sT(b), sB(b), sC(b), sD(b), eD(b));
    end
    fprintf('   (Ccoh at 1 kHz = %.2f)\n', exp(-w*1e6));
    semilogx(fBi, sB, '-o', 'LineWidth', 1.2, 'DisplayName', sprintf('coh. loss only, h_r = %g m', hr));
    semilogx(fBi, sD, '-s', 'LineWidth', 1.2, 'DisplayName', sprintf('both, h_r = %g m', hr));
end
semilogx(fBi, sT, 'k--', 'LineWidth', 1.4, 'DisplayName', 'target AF (h_r = 10 m)');
set(gca, 'XScale', 'log'); xlim([200 5000]);
xlabel('1/3-octave band (Hz)'); ylabel('std of level fluctuations (dB)'); legend('Location', 'northwest')
title('interaction with coherence loss')
