% VALIDATE_OSTASHEV_FIG3
% Reproduces Fig. 3 of Ostashev, Breton & Wilson, JASA 156(1), 107-118 (2024):
% coherence factor vs. frequency for three grazing angles, cloudy (QH = 0)
% and sunny (QH = 400 W/m^2) conditions.
%
% Setup of the paper (Sec. V B):
%   hr = 1.5 m, R1 = 200 m, hs = hr + R1*sin(alpha1),
%   alpha1 = 30 / 55 / 70 deg, Ts = 10 degC, zi = 1.5 km, ustar = 0.3 m/s.
%
% Expected behaviour: Ccoh decreases with frequency, increases with grazing
% angle (rho = 2*hr*cos(alpha2) gets smaller), and is lower for sunny than
% for cloudy conditions.

clear; close all; clc

hr     = 1.5;              % receiver height [m]
R1     = 200;              % direct path length [m]
alphas = [30 55 70];       % grazing angles [deg]
QHs    = [0 400];          % surface heat flux [W/m^2]

atm.ustar = 0.3;           % friction velocity [m/s]
atm.zi    = 1500;          % boundary layer height [m]
atm.Ts    = 10 + 273.15;   % surface temperature [K]

f = linspace(100, 10e3, 400);   % frequency [Hz]

figure('Position', [100 100 1000 400])
ls = {'-', '--', '-.'};

for iQ = 1:numel(QHs)
    atm.QH = QHs(iQ);
    subplot(1, 2, iQ); hold on; box on; grid on

    for ia = 1:numel(alphas)
        a1 = deg2rad(alphas(ia));
        hs = hr + R1*sin(a1);                       % source height [m]
        xr = R1*cos(a1);                            % horizontal range [m]
        a2 = atan((hs + hr)/xr);                    % Eq. (1)

        [w, Ccoh] = get_coherence_factor_slanted(hs, hr, a1, a2, atm, f);

        plot(f/1e3, Ccoh, ls{ia}, 'LineWidth', 1.4, ...
            'DisplayName', sprintf('\\alpha_1 = %d\\circ', alphas(ia)));

        fprintf('QH = %3d W/m^2 | a1 = %2d deg | hs = %6.1f m | a2 = %4.1f deg | w = %.3e\n', ...
            QHs(iQ), alphas(ia), hs, rad2deg(a2), w);
    end

    xlabel('Frequency (kHz)'); ylabel('Coherence factor')
    title(sprintf('Q_H = %d W/m^2', QHs(iQ)))
    ylim([0 1]); xlim([0 10]); legend('Location', 'northeast')
end

% ------------------------------------------------------------------
% Reference values from the Python prototype of the same equations
% (use these as a regression check of the MATLAB implementation):
%
%   QH =   0 W/m^2 : w = 9.344e-08 (30 deg), 4.159e-08 (55 deg), 1.787e-08 (70 deg)
%   QH = 400 W/m^2 : w = 3.687e-07 (30 deg), 1.762e-07 (55 deg), 7.582e-08 (70 deg)
%
% Corresponding coherence factors, e.g. QH = 400 W/m^2, alpha1 = 30 deg:
%   1 kHz -> 0.69, 2 kHz -> 0.23, 3 kHz -> 0.04
% and for QH = 0, alpha1 = 70 deg:
%   5 kHz -> 0.64, 7 kHz -> 0.42, 10 kHz -> 0.17
% ------------------------------------------------------------------