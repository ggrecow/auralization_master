function T = get_abl_turbulence(z, atm)
% GET_ABL_TURBULENCE  Variances and outer length scales of the von Karman
% turbulence spectrum at heights z, parameterized with Monin-Obukhov and
% mixed-layer similarity theory.
%
%   T = get_abl_turbulence(z, atm)
%
% Implements Eqs. (8)-(13) of
%   R. Pieren, D. Lincke, "Auralization of aircraft flyovers with
%   turbulence-induced coherence loss in ground effect", J. Acoust. Soc. Am.
%   151(4), 2453-2460 (2022),
% which follow Ostashev & Wilson, "Acoustics in Moving Inhomogeneous Media",
% 2nd ed. (2016). Shared by get_coherence_factor_slanted and
% get_coherence_factor_homogeneous, so that both use identical turbulence.
%
% INPUTS
%   z    : heights above ground [m] (any size; column vector recommended)
%   atm  : struct with
%            .ustar : friction velocity [m/s]
%            .QH    : surface sensible heat flux [W/m^2], positive upwards
%            .zi    : atmospheric boundary layer height [m]
%            .Ts    : surface temperature [K]
%
% OUTPUT
%   T    : struct with fields (same size as z where height-dependent)
%            .sigT2  : temperature variance [K^2]                     Eq. (8)
%            .LT     : temperature length scale [m]                   Eq. (11)
%            .sigvs2 : shear-produced velocity variance [m^2/s^2]     Eq. (9)
%            .Lvs    : shear-produced velocity length scale [m]       Eq. (12)
%            .sigvb2 : buoyancy-produced velocity variance [m^2/s^2]  Eq. (10)
%            .Lvb    : buoyancy-produced velocity length scale [m]    Eq. (13)
%
% LIMITATIONS
%   * Surface-layer quantities are valid for z < 0.1*zi; above that they are
%     frozen at their value at z = 0.1*zi (Pieren & Lincke, Sec. II B).
%   * The parameterization is for unstable (QH > 0) and neutral (QH = 0)
%     boundary layers. For stable conditions (QH < 0, typically at night)
%     Eqs. (8) and (11) become singular at z = Lo/10 and w* is undefined.
%     PROVISIONAL: QH < 0 is treated as neutral (QH = 0) and a warning is
%     issued. This keeps shear turbulence but drops temperature turbulence.
%
% Author: Gil Felix Greco, Braunschweig (27.10.2026)
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to 
%   develop this code under the author's direction. All codes were verified 
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent warned   % issue the stable-ABL warning only once per session

rho0 = 1.225;    % air density [kg/m^3]
cp   = 1005;     % specific heat at constant pressure [J/(kg K)]
grav = 9.81;     % gravitational acceleration [m/s^2]
kv   = 0.4;      % von Karman constant [-]

zc = min(z, 0.1*atm.zi);   % surface-layer validity limit
QH = atm.QH;

if QH < 0
    if isempty(warned)
        warning('coherenceLoss:stableBoundaryLayer', ...
            ['QH = %.1f W/m^2 < 0 (stable boundary layer) is outside the validity of the ' ...
             'similarity relations used here. It is treated as neutral (QH = 0).'], QH);
        warned = true;
    end
    QH = 0;
end

if QH > 0
    Lo    = -atm.ustar^3 * atm.Ts * rho0 * cp / (grav * kv * QH);   % Obukhov length (< 0)
    Tstar = QH / (rho0 * cp * atm.ustar);                           % temperature scale
    wstar = (atm.zi * grav * QH / (rho0 * cp * atm.Ts))^(1/3);      % convective velocity

    T.sigT2  = 4*Tstar^2 * (1 - 10*zc/Lo).^(-2/3);
    T.LT     = 2*zc .* (1 - 7*zc/Lo) ./ (1 - 10*zc/Lo);
    T.sigvb2 = 0.35 * wstar^2;
else
    T.sigT2  = zeros(size(zc));    % neutral: no temperature fluctuations
    T.LT     = 2*zc;
    T.sigvb2 = 0;                  % neutral: no buoyancy-produced turbulence
end

T.sigvs2 = 3 * atm.ustar^2;
T.Lvs    = 1.8 * zc;
T.Lvb    = 0.23 * atm.zi;

end
