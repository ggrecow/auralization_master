function [w, Ccoh] = get_coherence_factor_slanted(hs, hr, R1, alpha2, atm, f)
% GET_COHERENCE_FACTOR_SLANTED  Turbulence-induced coherence factor between
% the direct and the ground-reflected wave, for slanted propagation.
%
%   [w, Ccoh] = get_coherence_factor_slanted(hs, hr, R1, alpha2, atm, f)
%
% Implements Eq. (52) of
%   V. E. Ostashev, D. J. Breton, D. K. Wilson, "Interference of the direct
%   and ground-reflected waves in the atmosphere with volumetric
%   scattering", J. Acoust. Soc. Am. 156(1), 107-118 (2024),
%
%     Ccoh = exp{ -pi^2 k^2 R1 *
%                 Int_0^1 deta Int_0^inf Phi_eff(z(eta),0,kappa) *
%                 [1 - J0(rho*kappa*eta)] kappa dkappa }
%
% with z(eta) = hs - (hs-hr)*eta  (height along the direct path, source -> receiver)
% and  rho    = 2*hr*cos(alpha2)  (max. separation of the two paths, Eq. 48).
% The prefactor (hs-hr)/sin(alpha1) of Eq. (52) is written as the direct
% path length R1, to which it is identical. This avoids the singularity at
% alpha1 -> 0 and the sign problem for hs < hr.
%
% The exponent is proportional to k^2, so Ccoh(f) = exp(-w*f^2) with a
% single scalar w per source/receiver geometry.
%
% The von Karman spectrum (Eq. 53) is parameterized by get_abl_turbulence.
%
% INPUTS
%   hs      : source height above ground [m]
%   hr      : receiver height above ground [m]
%   R1      : direct (straight-line) source-receiver distance [m]
%   alpha2  : grazing angle of the reflected ray [rad], from the horizontal,
%             tan(alpha2) = (hs+hr)/xr (Ostashev Eq. 1)
%   atm     : struct with .ustar [m/s], .QH [W/m^2], .zi [m], .Ts [K],
%             and optionally .c0 [m/s] (default 20.05*sqrt(Ts))
%   f       : (optional) frequency vector [Hz] for the second output
%
% OUTPUTS
%   w       : scalar, Ccoh = exp(-w*f^2)   [s^2]
%   Ccoh    : coherence factor at f (empty if f not given)
%
% VALIDITY
%   * Derived for slanted propagation: hs >> hr and hs >~ xr/2, i.e. grazing
%     angles >~ 25 deg. For small angles and hs >> hr the formula
%     approximately coincides with the near-horizontal quasi-homogeneous
%     result (Ostashev Eq. 46), but it is not derived for that regime, and
%     it is not equivalent to the homogeneous model of Pieren & Lincke.
%   * Neglects the leg from the reflection point to the receiver; requires a
%     receiver within a few metres of the ground.
%   * hr = 0 gives rho = 0 and w = 0 (no second path).
%
% VALIDATION
%   Reproduces Fig. 3 of Ostashev et al. (2024), see validate_ostashev.m.
%
% Author: Gil Felix Greco, Braunschweig (27.10.2026)
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to 
%   develop this code under the author's direction. All codes were verified 
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(5, 6);
if nargin < 6, f = []; end

if ~isfield(atm, 'c0') || isempty(atm.c0)
    atm.c0 = 20.05*sqrt(atm.Ts);
end
c0 = atm.c0;

% ---------------- geometry ----------------
rho = 2*hr*cos(alpha2);                    % Eq. (48): path separation [m]

if hr <= 0 || R1 <= 0 || rho <= 0
    w = 0;
    Ccoh = ones(size(f));
    if isempty(f), Ccoh = []; end
    return
end

% ---------------- integration grids ----------------
nEta = 200;                                % along the propagation path
nKap = 400;                                % turbulence wavenumber
eta  = linspace(0, 1, nEta).';             % column
kap  = logspace(-4, 2, nKap);              % row, [rad/m]
lkap = log(kap);

z = hs - (hs - hr)*eta;                    % height along the path [m]

% ---------------- effective von Karman spectrum, Eq. (53) ----------------
T = get_abl_turbulence(z, atm);

A  = gamma(11/6) / (pi^1.5 * gamma(1/3));
t1 = T.sigT2 .* T.LT.^3 ./ (atm.Ts^2 * (1 + kap.^2 .* T.LT.^2).^(11/6));
t2 = (22/3) * T.sigvs2 * T.Lvs.^5 .* kap.^2 ./ (c0^2 * (1 + kap.^2 .* T.Lvs.^2).^(17/6));
t3 = (22/3) * T.sigvb2 * T.Lvb^5  .* kap.^2 ./ (c0^2 * (1 + kap.^2 *  T.Lvb^2 ).^(17/6));
Phi = A * (t1 + t2 + t3);                  % [nEta x nKap]

% ---------------- double integral ----------------
integrand = Phi .* (1 - besselj(0, rho*(kap.*eta))) .* kap;
Ikap      = trapz(lkap, integrand .* kap, 2);   % dkappa = kappa dln(kappa)
I         = trapz(eta, Ikap);

% ---------------- exponent, with k^2 = (2*pi*f/c0)^2 pulled out ----------------
w = 4*pi^4 * R1 * I / c0^2;

if ~isempty(f)
    Ccoh = reshape(exp(-w * f(:).^2), size(f));
else
    Ccoh = [];
end

end
