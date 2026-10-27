% COMPARE_COHERENCE_MODELS
% Compares the two coherence-factor formulations along a flyover trajectory:
%
%   (a) Pieren & Lincke (2022), Eq. 4  -> get_coherence_factor_homogeneous
%       homogeneous turbulence, all parameters evaluated at z = h ~ 2*hr,
%       applied along the whole path.
%
%   (b) Ostashev et al. (2024), Eq. 52 -> get_coherence_factor_slanted
%       quasi-homogeneous turbulence, spectrum integrated from hr up to hs,
%       separation rho = 2*hr*cos(alpha2).
%
% Since the turbulence length scales grow with height, (a) is expected to
% overestimate the decorrelation when the aircraft is high above the
% receiver, i.e. near the closest point of approach. This script quantifies
% that difference.
%
% Replace the synthetic trajectory below by your own PANAM trajectory and
% receiver position to get the numbers relevant for your cases.

clear; close all; clc

% ---------------- trajectory (straight and level flyover) ----------------
hs_const = 300;                       % aircraft altitude [m]
hr       = 1.5;                       % receiver height [m]
x        = linspace(-4000, 4000, 401).';   % along-track position [m]

% ---------------- atmosphere ----------------
atm.ustar = 0.3;                      % friction velocity [m/s]
atm.QH    = 150;                      % surface heat flux [W/m^2]
atm.zi    = 1000;                     % boundary layer height [m]
atm.Ts    = 15 + 273.15;              % surface temperature [K]

f = [250 500 1000 2000].';            % frequencies to display [Hz]

% ---------------- loop over trajectory ----------------
n   = numel(x);
wH  = zeros(n,1);   % homogeneous
wS  = zeros(n,1);   % slanted
a1  = zeros(n,1);

for i = 1:n
    hs  = hs_const;
    xr  = abs(x(i));                             % horizontal range [m]
    R1  = sqrt(xr^2 + (hs-hr)^2);                % direct path length [m]

    alpha1 = asin( (hs-hr)/R1 );                 % grazing angle, direct
    alpha2 = atan( (hs+hr)/max(xr, eps) );       % grazing angle, reflected (Eq. 1)
    % (alpha1 only used for plotting; no guard needed with the R1 prefactor)

    wH(i) = get_coherence_factor_homogeneous(hs, hr, R1, atm);
    wS(i) = get_coherence_factor_slanted(hs, hr, R1, alpha2, atm);
    a1(i) = rad2deg(alpha1);
end

% ---------------- plots ----------------
figure('Position', [80 80 1200 700])

subplot(2,2,1)
semilogy(x/1e3, wH, 'LineWidth', 1.3); hold on; grid on; box on
semilogy(x/1e3, wS, 'LineWidth', 1.3)
xlabel('Along-track position (km)'); ylabel('w (s^2)')
legend('homogeneous (Pieren & Lincke)', 'slanted (Ostashev)', 'Location', 'south')
title('coherence parameter, C_{coh} = exp(-w f^2)')

subplot(2,2,2)
plot(x/1e3, wH./wS, 'LineWidth', 1.3); grid on; box on
xlabel('Along-track position (km)'); ylabel('w_{hom} / w_{slant}')
title('ratio of the two formulations')

subplot(2,2,3); hold on; grid on; box on
co = lines(numel(f));
for k = 1:numel(f)
    plot(x/1e3, exp(-wH*f(k)^2), '-',  'Color', co(k,:), 'LineWidth', 1.3, ...
        'DisplayName', sprintf('%d Hz, hom.', f(k)))
    plot(x/1e3, exp(-wS*f(k)^2), '--', 'Color', co(k,:), 'LineWidth', 1.3, ...
        'DisplayName', sprintf('%d Hz, slant.', f(k)))
end
xlabel('Along-track position (km)'); ylabel('C_{coh}')
ylim([0 1]); legend('Location', 'south', 'NumColumns', 2)
title('coherence factor over the flyover')

subplot(2,2,4)
plot(x/1e3, a1, 'LineWidth', 1.3); grid on; box on
xlabel('Along-track position (km)'); ylabel('\alpha_1 (deg)')
title('grazing angle of the direct ray')

% ---------------- summary ----------------
[~, iOver] = min(abs(x));
fprintf('Overhead (alpha1 = %.0f deg): w_hom = %.3e, w_slant = %.3e, ratio = %.2f\n', ...
    a1(iOver), wH(iOver), wS(iOver), wH(iOver)/wS(iOver));
fprintf('Far field (alpha1 = %.0f deg): w_hom = %.3e, w_slant = %.3e, ratio = %.2f\n', ...
    a1(1), wH(1), wS(1), wH(1)/wS(1));
fprintf('\nAt 1 kHz overhead:  C_hom = %.3f, C_slant = %.3f\n', ...
    exp(-wH(iOver)*1e6), exp(-wS(iOver)*1e6));
