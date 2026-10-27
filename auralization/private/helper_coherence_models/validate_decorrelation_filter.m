% VALIDATE_DECORRELATION_FILTER
% Checks that get_decorrelation_filter reproduces its target coherence factor
% in the ensemble mean, i.e. |E[B(f)]| = Ccoh(f), and inspects the magnitude
% response and the impulse response (cf. Fig. 2 of Pieren & Lincke 2022).

clear; close all; clc

fs   = 48e3;
nfft = 2^14;                  % as used in MASTER_auralization_EngineAirframe
nMC  = 400;                   % Monte-Carlo realizations
rs   = RandStream('mt19937ar', 'Seed', 1);

% three representative values from the Ostashev Fig. 3 geometry
wList = [3.687e-7, 9.344e-8, 1.787e-8];
lab   = {'sunny, \alpha_1 = 30\circ', 'cloudy, \alpha_1 = 30\circ', 'cloudy, \alpha_1 = 70\circ'};

nBins = nfft/2 + 1;
f     = (0:nBins-1).' * fs/nfft;

figure('Position', [100 100 1100 400])

for iw = 1:numel(wList)
    w = wList(iw);

    Bsum  = zeros(nBins, 1);
    Bmag  = zeros(nBins, 1);
    for m = 1:nMC
        [B, bhat, Ccoh] = get_decorrelation_filter(w, nfft, fs, rs);
        Bsum = Bsum + B;
        Bmag = Bmag + abs(B);
    end
    Cest = abs(Bsum/nMC);      % achieved coherence
    Bmag = Bmag/nMC;           % mean magnitude response

    subplot(1,3,1); hold on; box on; grid on
    plot(f/1e3, Ccoh, 'k--', 'LineWidth', 1.2)
    plot(f/1e3, Cest, 'LineWidth', 1.2, 'DisplayName', lab{iw})
    xlim([0 8]); ylim([0 1])
    xlabel('Frequency (kHz)'); ylabel('Coherence factor')
    title('target (dashed) vs. achieved')

    subplot(1,3,2); hold on; box on; grid on
    plot(f/1e3, 20*log10(Bmag), 'LineWidth', 1.2)
    xlim([0 8]); ylim([-3 1])
    xlabel('Frequency (kHz)'); ylabel('mean |B| (dB)')
    title('magnitude response (all-pass check)')

    % report at a few frequencies
    for ff = [500 1000 2000 4000]
        i = round(ff/fs*nfft) + 1;
        fprintf('w = %.3e | %5d Hz | target %.3f | achieved %.3f | |B| %.3f\n', ...
            w, ff, Ccoh(i), Cest(i), Bmag(i));
    end
    fprintf('\n');
end

% impulse response of the strongest decorrelation case, cf. their Fig. 2(b)
[~, bhat] = get_decorrelation_filter(wList(1), nfft, fs, rs);
subplot(1,3,3)
t = ((0:numel(bhat)-1).' - numel(bhat)/2)/fs*1e3;
plot(t, bhat, 'LineWidth', 0.8); grid on; box on
xlabel('Time (ms)'); ylabel('$\hat{b}$', 'Interpreter', 'latex')
title('impulse response')

% ------------------------------------------------------------------
% Reference from the Python prototype of the same design (nMC = 400):
%   w = 3.7e-07 : 1 kHz target 0.691 -> 0.703 | 2 kHz 0.227 -> 0.200
%   w = 9.3e-08 : 1 kHz target 0.911 -> 0.932 | 2 kHz 0.689 -> 0.686
%   w = 1.8e-08 : 2 kHz target 0.930 -> 0.952 | 4 kHz 0.750 -> 0.769
% Mean |B| stays within about -1.4 dB, reaching that only where Ccoh -> 0.
% ------------------------------------------------------------------
