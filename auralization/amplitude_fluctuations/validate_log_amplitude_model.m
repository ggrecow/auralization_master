% VALIDATE_LOG_AMPLITUDE_MODEL
% Checks the semi-empirical model for turbulence-induced level fluctuations
% (get_log_amplitude_fluctuation) against the values published by Lincke &
% Pieren (2024) [LP24] and Lincke, Kawai & Pieren (2025) [LKP25]:
%
%   1) model curves sigma(f) of LKP25 Fig. 2 (cases I, III, IV, VI: A320 at
%      hs = 700 m and A340 at hs = 400 m, receiver at hr = 10 m, overhead,
%      meteo conditions A and C of LKP25 Table 1), for three variants:
%        - DEFAULT: saturation of the total, effective distance Leff (LP24 Eq. 41)
%        - saturation of the total, full path length (no Leff)
%        - saturation per component (LP24 Eq. 37 as printed), Leff
%   2) limit frequency f_lim for the six cases of LKP25 (Table 2), for information
%   3) level of the fluctuations for the two example conditions of LP24
%      (Table 2, Sec. 7): 1.1 dB (condition I) and 2.7 dB (condition II)
%   4) dependence on distance for a level flyover, per 1/3-octave band,
%      compared with the range measured at Zurich airport (LP24 Sec. 4.1:
%      sigma = 1...4 dB, mean 2.2 dB, 446-1120 Hz, direct distances 570-3250 m)
%
% The LKP25 Fig. 2 values are read off the figure (approx. +-0.1 dB).
%
% EXPECTED OUTCOME (18.11.2026, GNU Octave 8.4)
%   1) rms difference to LKP25 Fig. 2 (I / III / IV / VI):
%        DEFAULT (total, Leff)        0.98 / 0.23 / 0.78 / 0.09 dB
%        total, full path             0.07 / 0.23 / 0.14 / 0.09 dB
%        per component, Leff          0.98 / 0.93 / 0.78 / 0.59 dB
%      Condition C (III, VI; zi > hs, Leff not involved) is reproduced with the
%      total weighting. Condition A (I, IV; zi = 52 m < hs) is reproduced only
%      without Leff, but see 3): LP24 condition I agrees only WITH Leff, so the
%      papers disagree here; Leff (LP24, fitted to measurements) is the DEFAULT.
%      Condition A has QH < 0 (stable), treated as neutral by get_abl_turbulence
%      (warning); this changes the result by less than 0.2 dB.
%   2) f_lim here (DEFAULT): I and IV NaN (> 10 kHz), II 4057, III 416,
%      V 7327, VI 1013 Hz; LKP25: 1373, 620, 155, 3039, 1370, 230 Hz. Not
%      reproduced, although sigma(f) of III and VI matches Fig. 2, i.e. Table 2
%      was obtained differently; not used as a reference.
%   3) LP24 (446-1120 Hz, overhead): I 0.8 dB (full path 2.2 dB), II 3.0 dB.
%      The paper values (1.1 / 2.7 dB) are averages over +-10 s around the
%      closest point; the flight path of the examples is not given, so only
%      the order can be compared.
%   4) strongly convective conditions (LP24 II, LKP25 C): 446-1120 Hz value
%      from ~2.9-3.0 dB (overhead) to ~3.5 dB at 3 km, i.e. within the measured
%      range of 1-4 dB. The total weighting saturates at <chi^2> = b^2 = 1/6,
%      i.e. 3.5 dB. The weak-scattering theory alone (dashed) is much higher,
%      as reported in LP24 (Sec. 6.2, mean absolute error 2.9 dB).
%
% Author: Gil Felix Greco, Braunschweig (29.10.2026)
% 18.11.2026: comparison with LKP25 Fig. 2 (saturation of the total, effective distance)
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to
%   develop this code under the author's direction. All codes were verified
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

hr = 10;

%            u* [m/s]  QH [W/m2]  zi [m]  Ts [degC]   (LKP25 Table 1)
meteo = struct( ...
    'name',  {'A (low)', 'B (intermediate)', 'C (high)'}, ...
    'ustar', {0.10,  0.11, 0.28}, ...
    'QH',    {-7,    22,   109}, ...
    'zi',    {52,    100,  1040}, ...
    'Ts',    {16.3,  19.0, 22.0} );

%% 1) LKP25 Fig. 2: model curves sigma(f), cases I, III, IV, VI

fRef = [50 100 200 500 1000 2000 5000 10000];
fig2 = struct( ...                                % read off LKP25 Fig. 2 (approx. +-0.1 dB)
    'id',    {'I', 'III', 'IV', 'VI'}, ...
    'hs',    {700, 700, 400, 400}, ...
    'meteo', {1, 3, 1, 3}, ...
    'sigma', {[0.28 0.50 0.85 1.50 2.20 2.70 3.10 3.25], ...
              [1.00 1.35 1.80 2.30 2.60 3.00 3.35 3.45], ...
              [0.20 0.35 0.55 1.05 1.50 2.05 2.60 2.80], ...
              [0.62 0.95 1.45 2.15 2.50 2.90 3.30 3.42]} );

variants = { 'DEFAULT (total, L_{eff})', struct(); ...
             'total, full path',         struct('effectiveDistance', false); ...
             'per component, L_{eff}',   struct('saturationMode', 'component') };
styles = {'-', '--', ':'};

f = logspace(log10(50), log10(10e3), 200);

figure('Position', [100 100 1200 700])
fprintf('1) rms difference to LKP25 Fig. 2 [dB]\n');
fprintf('   %-26s |  I    |  III  |  IV   |  VI\n', 'variant');
for v = 1:size(variants, 1)
    rmsd = zeros(1, numel(fig2));
    for c = 1:numel(fig2)
        m   = meteo(fig2(c).meteo);
        atm = struct('ustar', m.ustar, 'QH', m.QH, 'zi', m.zi, 'Ts', m.Ts + 273.15);
        hs  = fig2(c).hs;
        sRef = get_log_amplitude_fluctuation(hs, hr, hs - hr, atm, fRef, variants{v,2});
        rmsd(c) = sqrt(mean((sRef - fig2(c).sigma).^2));

        sig = get_log_amplitude_fluctuation(hs, hr, hs - hr, atm, f, variants{v,2});
        subplot(2, 2, c); hold on; box on; grid on
        semilogx(f, sig, styles{v}, 'LineWidth', 1.4, 'DisplayName', variants{v,1});
        if v == 1
            semilogx(fRef, fig2(c).sigma, 'ko', 'LineWidth', 1.2, 'DisplayName', 'LKP25 Fig. 2');
        end
        set(gca, 'XScale', 'log'); xlim([50 10e3]); ylim([0 5])
        xlabel('Frequency (Hz)'); ylabel('\sigma_\chi (dB)')
        title(sprintf('case %s: h_s = %d m, meteo %s', fig2(c).id, hs, m.name))
    end
    fprintf('   %-26s | %5.2f | %5.2f | %5.2f | %5.2f\n', strrep(variants{v,1}, '_{eff}', '_eff'), rmsd);
end
subplot(2, 2, 1); legend('Location', 'northwest')

%% 2) LKP25 Table 2: limit frequency (information only)

cases = struct( ...
    'id',    {'I', 'II', 'III', 'IV', 'V', 'VI'}, ...
    'hs',    {700, 700, 700, 400, 400, 400}, ...
    'meteo', {1, 2, 3, 1, 2, 3}, ...
    'flim_paper', {1373, 620, 155, 3039, 1370, 230} );   % LKP25 Table 2 [Hz]

fFine = logspace(log10(20), log10(10e3), 1000);
fprintf('\n2) limit frequency (DEFAULT), overhead, hr = %g m\n', hr);
fprintf('   case | aircraft hs | meteo            | f_lim here | f_lim LKP25 | sigma at 1 kHz\n');
for c = 1:numel(cases)
    m   = meteo(cases(c).meteo);
    atm = struct('ustar', m.ustar, 'QH', m.QH, 'zi', m.zi, 'Ts', m.Ts + 273.15);
    hs  = cases(c).hs;
    [~, out] = get_log_amplitude_fluctuation(hs, hr, hs - hr, atm, fFine);
    s1k = get_log_amplitude_fluctuation(hs, hr, hs - hr, atm, 1000);
    fprintf('   %4s | %11d | %-16s | %8.0f Hz | %8.0f Hz  | %.1f dB\n', ...
        cases(c).id, hs, m.name, out.flim, cases(c).flim_paper, s1k);
end

%% 3) LP24: example conditions I and II (A320 at ~700 m, hr = 10 m)

%                  u*     QH      zi     Ts [degC]     sigma LP24 (+-10 s average)
condLP24 = { 'I',  0.09,  10.6,   48.5,  288.3-273.15, 1.1; ...
             'II', 0.12,  275.9,  727,   295.8-273.15, 2.7 };

fBand = linspace(446, 1120, 50);   % LP24 Sec. 5.2: average over linearly spaced frequencies
fprintf('\n3) LP24 example conditions (overhead, hs = 700 m, hr = 10 m, 446-1120 Hz)\n');
for c = 1:size(condLP24, 1)
    atm = struct('ustar', condLP24{c,2}, 'QH', condLP24{c,3}, 'zi', condLP24{c,4}, 'Ts', condLP24{c,5} + 273.15);
    s   = get_log_amplitude_fluctuation(700, 10, 690, atm, fBand);
    sF  = get_log_amplitude_fluctuation(700, 10, 690, atm, fBand, struct('effectiveDistance', false));
    sW  = get_log_amplitude_fluctuation(700, 10, 690, atm, fBand, struct('saturation', false));
    fprintf('   condition %-2s: DEFAULT %.1f dB | full path %.1f dB | weak scattering %.1f dB | LP24 %.1f dB (+-10 s average)\n', ...
        condLP24{c,1}, sqrt(mean(s.^2)), sqrt(mean(sF.^2)), sqrt(mean(sW.^2)), condLP24{c,6});
end

%% 4) distance dependence, level flyover at hs = 700 m, hr = 10 m (DEFAULT)

xr    = linspace(0, 3200, 120);         % horizontal distance [m]
hs    = 700;
R1    = sqrt(xr.^2 + (hs - hr)^2);
fOB   = [250 500 1000 2000 4000];       % 1/3-octave band centres [Hz]

figure('Position', [100 100 1100 420])
fprintf('\n4) 446-1120 Hz vs. direct distance (DEFAULT)\n');
for c = 1:2
    if c == 1
        atm = struct('ustar', 0.12, 'QH', 275.9, 'zi', 727, 'Ts', 295.8); ttl = 'LP24 condition II (zi = 727 m)';
    else
        atm = struct('ustar', 0.28, 'QH', 109, 'zi', 1040, 'Ts', 295.15); ttl = 'LKP25 condition C (zi = 1040 m)';
    end
    S  = zeros(numel(R1), numel(fOB));
    SB = zeros(numel(R1), 1); SW = SB;
    for i = 1:numel(R1)
        S(i,:) = get_log_amplitude_fluctuation(hs, hr, R1(i), atm, fOB);
        SB(i)  = sqrt(mean(get_log_amplitude_fluctuation(hs, hr, R1(i), atm, fBand).^2));
        SW(i)  = sqrt(mean(get_log_amplitude_fluctuation(hs, hr, R1(i), atm, fBand, struct('saturation', false)).^2));
    end
    fprintf('   %s: %.1f dB at %.0f m, %.1f dB at %.0f m\n', ttl, SB(1), R1(1), SB(end), R1(end));

    subplot(1, 2, c); hold on; box on; grid on
    patch([570 3250 3250 570], [1 1 4 4], [0.9 0.9 0.9], 'EdgeColor', 'none', ...
        'DisplayName', 'LP24 measured range (446-1120 Hz)');
    for k = 1:numel(fOB)
        plot(R1, S(:,k), 'LineWidth', 1.2, 'DisplayName', sprintf('%g Hz', fOB(k)));
    end
    plot(R1, SB, 'k-', 'LineWidth', 1.8, 'DisplayName', '446-1120 Hz');
    plot(R1, SW, 'k--', 'LineWidth', 1.0, 'DisplayName', '446-1120 Hz, weak scattering');
    xlabel('Direct distance R_1 (m)'); ylabel('\sigma_\chi (dB)'); ylim([0 8]); xlim([R1(1) R1(end)])
    title(ttl); legend('Location', 'northwest')
end
