function [w, Ccoh] = get_coherence_factor_homogeneous(hs, hr, r, atm, f)
% GET_COHERENCE_FACTOR_HOMOGENEOUS  Turbulence-induced coherence factor
% between direct and ground-reflected wave, homogeneous formulation.
%
%   [w, Ccoh] = get_coherence_factor_homogeneous(hs, hr, r, atm, f)
%
% Implements Eq. (4) of
%   R. Pieren, D. Lincke, "Auralization of aircraft flyovers with
%   turbulence-induced coherence loss in ground effect", J. Acoust. Soc. Am.
%   151(4), 2453-2460 (2022),
% i.e. the closed-form solution of Ostashev & Wilson (2016) for statistically
% homogeneous von Karman turbulence and near-horizontal propagation:
%
%   Ccoh = exp{ -SUM_i (2 gamma_i r L_i / h) * Int_0^{h/L_i} A_i(eta) deta }
%
% with h = 2/(1/hr + 1/hs) (Eq. 14), kernels A (Eq. 5, temperature) and
% A' (Eq. 6, velocity), and extinction coefficients gamma_i (Eq. 7).
%
% All turbulence quantities are evaluated at the single height z = h and
% applied along the whole path of length r. Use this function to reproduce
% Pieren & Lincke (2022), and as reference for get_coherence_factor_slanted.
%
% INPUTS
%   hs, hr : source and receiver heights above ground [m]
%   r      : source-receiver distance [m]
%   atm    : struct with .ustar, .QH, .zi, .Ts, optionally .c0
%   f      : (optional) frequency vector [Hz]
%
% OUTPUTS
%   w      : scalar, Ccoh = exp(-w*f^2) [s^2]
%   Ccoh   : coherence factor at f (empty if f not given)
%
% VALIDITY
%   Near-horizontal propagation (grazing angles <~ 20 deg), homogeneous
%   turbulence. For elevated sources it applies near-ground turbulence along
%   the whole path and overestimates the decorrelation.
%
% Author: Gil Felix Greco, Braunschweig (27.10.2026)
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to 
%   develop this code under the author's direction. All codes were verified 
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(4, 5);
if nargin < 5, f = []; end

if ~isfield(atm, 'c0') || isempty(atm.c0)
    atm.c0 = 20.05*sqrt(atm.Ts);
end
c0 = atm.c0;

% ---- maximum transverse separation, Eq. (14) ----
if hr <= 0 || hs <= 0 || r <= 0
    w = 0;
    Ccoh = ones(size(f));
    if isempty(f), Ccoh = []; end
    return
end
h = 2 / (1/hr + 1/hs);

% ---- turbulence at z = h, Eqs. (8)-(13) ----
T = get_abl_turbulence(h, atm);

% ---- extinction coefficients per k^2, Eq. (7) ----
pre    = sqrt(pi) * gamma(5/6) / (4 * gamma(1/3));
gamT   = pre *     T.sigT2  * T.LT  / atm.Ts^2;
gamvs  = pre * 4 * T.sigvs2 * T.Lvs / c0^2;
gamvb  = pre * 4 * T.sigvb2 * T.Lvb / c0^2;

% ---- the three path integrals, Eq. (4) ----
expo_k2 =   il_segment(gamT,  T.LT,  h, r, @il_A) ...
          + il_segment(gamvs, T.Lvs, h, r, @il_Aprime) ...
          + il_segment(gamvb, T.Lvb, h, r, @il_Aprime);

% ---- pull k^2 = (2*pi*f/c0)^2 out ----
w = expo_k2 * (2*pi/c0)^2;

if ~isempty(f)
    Ccoh = reshape(exp(-w * f(:).^2), size(f));
else
    Ccoh = [];
end

end

% =====================================================================
function S = il_segment(gam, L, h, r, kernel)
% 2*gamma*r*L/h * Int_0^{h/L} kernel(eta) deta
if gam == 0 || L == 0
    S = 0; return
end
eta = linspace(1e-9, h/L, 4000);
S   = 2 * gam * r * L / h * trapz(eta, kernel(eta));
end

% =====================================================================
function y = il_A(eta)      % Eq. (5), temperature term
y = 1 - 2^(1/6) * eta.^(5/6) ./ gamma(5/6) .* besselk(5/6, eta);
end

% =====================================================================
function y = il_Aprime(eta) % Eq. (6), velocity terms
y = 1 - 2^(1/6) * eta.^(5/6) ./ gamma(5/6) .* ...
    ( besselk(5/6, eta) - 0.5*eta .* besselk(1/6, eta) );
end
