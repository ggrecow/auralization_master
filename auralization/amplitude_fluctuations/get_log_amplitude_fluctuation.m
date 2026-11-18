function [sigma_dB, out] = get_log_amplitude_fluctuation(hs, hr, R1, atm, f, opts)
% GET_LOG_AMPLITUDE_FLUCTUATION  Standard deviation of turbulence-induced
% level (log-amplitude) fluctuations for slanted propagation from an elevated
% source, following the semi-empirical model of Lincke & Pieren (2024).
%
%   [sigma_dB, out] = get_log_amplitude_fluctuation(hs, hr, R1, atm, f, opts)
%
% The model takes the weak-scattering (Rytov) theory of Ostashev & Wilson,
% written in terms of strength parameters Phi^2 and wave parameters D, and
% adds an empirical saturation weighting fitted to > 5000 aircraft flyovers:
%
%   <chi^2>_i  = Phi_i^2/2 * (1 - Hbar_i(D_i))                 weak scattering, LP24 Eq. (28)
%   Phi_max^2  = 2 b^2 / sum_i beta_i (1 - Hbar_i(D_i))        LP24 Eq. (30), b^2 = 1/6 (Eq. 43)
%   <chi^2>^   = sum_i <chi^2>_i / (1 + max_i(Phi_i^2)/Phi_max^2)   saturation (see below)
%   sigma_dB   = 20*log10(e) * sqrt(<chi^2>^)                  LP24 Eq. (3)
%
% with i = T (temperature), vs (shear-produced velocity), vb (buoyancy-
% produced velocity), and the strength and wave parameters integrated along
% the direct path (LP24 Eqs. 34-36, LKP25 Eqs. 3, 4 and 8):
%
%   Phi_T^2 = sqrt(pi)G(a-1)/(2G(a-3/2)) k^2 L Int sigT^2 LT / T0^2 deta
%   Phi_v^2 = 2 sqrt(pi)G(a-1)/G(a-3/2)  k^2 L Int sigv^2 Lv / c0^2 deta
%   D_i     = L/k Int 1/L_i^2 deta
%
% For sources above the boundary layer (hs > zi), L is replaced by the
% effective propagation distance (LP24 Eq. 41, q = 0.6)
%
%   Leff = L_ABL * (zi/hs)^(-q(1 - zi/hs))
%
% where L_ABL is the part of the direct path below zi.
%
% The turbulence (variances and length scales) comes from get_abl_turbulence,
% i.e. it is IDENTICAL to the one used for the coherence loss.
%
% INPUTS
%   hs    : source height above ground [m]
%   hr    : receiver height above ground [m]
%   R1    : direct (straight-line) source-receiver distance [m]
%   atm   : struct with .ustar [m/s], .QH [W/m^2], .zi [m], .Ts [K],
%           optionally .c0 [m/s] (default 20.05*sqrt(Ts)). Same struct as for
%           get_coherence_factor_slanted.
%   f     : frequency vector [Hz]
%   opts  : (optional) struct
%       .b2         weak/strong scattering boundary, default 1/6 (LP24 Eq. 43)
%       .q          effective-distance exponent, default 0.6 (LP24 Eq. 41)
%       .saturation true (default): semi-empirical model
%                   false: weak-scattering theory only, Eq. (28)
%       .saturationMode  'total' (DEFAULT): weighting of the summed <chi^2>
%                   with the strongest strength parameter (reproduces LKP25 Fig. 2)
%                   'component': per turbulence component, LP24 Eq. (37) as printed
%       .effectiveDistance true (DEFAULT): Leff for hs > zi (LP24 Eq. 41);
%                   false: full path length (as in LKP25 Fig. 2, see below)
%
% OUTPUTS
%   sigma_dB : standard deviation of the level fluctuations at f [dB], same size as f
%   out      : struct with intermediate results (vectors over f, columns T, vs, vb)
%       .phi2      strength parameters Phi^2           [nf x 3]
%       .D         path-averaged wave parameters D (LP24 Eq. 34, info only) [nf x 3]
%       .omH       effective 1 - Hbar, see below       [nf x 3]
%       .phi2max   Phi_max^2                           [nf x 1]
%       .chi2weak  weak-scattering <chi^2>, Eq. (28)   [nf x 1]
%       .chi2      semi-empirical <chi^2>^             [nf x 1]
%       .flim      lowest f with max(Phi^2) >= Phi_max^2 (LKP25 Eq. 10), NaN if none
%       .Leff      (effective) propagation distance [m]
%
% IMPLEMENTATION CHOICES (not unambiguously defined in the papers)
%   * <chi^2>_i follows LP24 Eq. (28) literally: the diffraction term
%     1 - H(D(z), eta) is evaluated with the LOCAL wave parameter
%     D(z) = L/(k L_i(z)^2) inside the path integral. The resulting
%     effective 1 - Hbar_i = Int Phi_i^2-density (1 - H) deta / Phi_i^2 is
%     also used in Phi_max^2 (LP24 Eq. 30, written there with H(D) of a
%     single, path-integrated D). Using a path-averaged D instead (LP24
%     Eq. 34) makes even the weak-scattering result depend strongly on the
%     receiver height (D is dominated by the small eddies near the ground):
%     at 1 kHz, condition C of LKP25, 12.5 dB at hr = 1.2 m vs 6.6 dB at
%     hr = 10 m, which is unphysical. The local form gives 6.0 dB for both.
%   * Saturation: LP24 Eq. (37) and LKP25 Eq. (2) print the weighting per
%     turbulence component, 1/(1 + Phi_i^2/Phi_max^2) (LP24 Eq. (42) even with
%     Phi/Phi_max). Applied per component, the sum can exceed b^2 by far
%     (up to ~5 dB at high frequencies and long distances, above the 1-4 dB
%     measured in LP24). Applied to the TOTAL, with the strongest strength
%     parameter, the result saturates at <chi^2> = b^2 (3.5 dB) and reproduces
%     the model curves published in LKP25 Fig. 2 for condition C (rms
%     difference 0.09 dB for case VI, 0.23 dB for case III, vs. 0.59 and
%     0.93 dB per component). The total weighting is therefore the DEFAULT.
%   * Effective distance (hs > zi): the papers disagree. LKP25 Fig. 2 (cases I
%     and IV, zi = 52 m) is reproduced only with the full path (rms 0.07 and
%     0.14 dB; with Leff 0.98 and 0.78 dB), whereas the example condition I of
%     LP24 (zi = 48.5 m, 1.1 dB) agrees only with Leff (0.85 dB at overhead,
%     full path 2.2 dB). Leff is kept as DEFAULT, since LP24 is the model
%     fitted to the measurements (mean absolute error 0.4 dB).
%   * The limit frequencies f_lim of LKP25 Table 2 are not reproduced (e.g.
%     case III: ~420 Hz here vs. 155 Hz), although sigma(f) of case III
%     matches LKP25 Fig. 2; they are not used as a reference.
%     See validate_log_amplitude_model.m. All choices are made in this
%     function only, so they can be changed once clarified.
%   * The corrections of LP24 for turbulence decay in the afternoon (Eq. 38)
%     and low-level jets (Eq. 40) need inputs that are not used here. For the
%     former, the user may enter QH(t - 3t*) as surface heat flux.
%   * hr <= 0 (flush-mounted receiver): D(z) -> Inf at the ground, where
%     the limit 1 - H = 1 is used.
%
% VALIDITY (LP24 Sec. 6.5)
%   Fitted to narrow-body airliner flyovers within +-10 s of the closest point,
%   direct distances 570-3250 m, 446-1120 Hz, no rain, T > 0 degC, wind < 7 m/s,
%   ERA5 turbulence parameters (do not use in-situ u*, LP24 Fig. 7).
%   Turbulence parameterization valid for QH >= 0 only (see get_abl_turbulence).
%
% References
%   [LP24]  D. Lincke, R. Pieren, Acta Acustica 8, 47 (2024). doi:10.1051/aacus/2024036
%   [LKP25] D. Lincke, C. Kawai, R. Pieren, Acta Acustica 9, 69 (2025). doi:10.1051/aacus/2025053
%   V. E. Ostashev, D. K. Wilson, Acta Acust. united Ac. 104, 571-585 (2018)
%   V. E. Ostashev, D. K. Wilson, JASA 141(3), 2079-2092 (2017)
%
% Author: Gil Felix Greco, Braunschweig (29.10.2026)
% 18.11.2026: saturation applied to the total (DEFAULT, reproduces LKP25 Fig. 2),
% options saturationMode and effectiveDistance
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to 
%   develop this code under the author's direction. All codes were verified 
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(5, 6);
if nargin < 6, opts = struct; end
if ~isfield(opts, 'b2'),         opts.b2 = 1/6;         end
if ~isfield(opts, 'q'),          opts.q = 0.6;          end
if ~isfield(opts, 'saturation'), opts.saturation = true; end
if ~isfield(opts, 'saturationMode'),    opts.saturationMode = 'total';  end
if ~isfield(opts, 'effectiveDistance'), opts.effectiveDistance = true;  end

if ~isfield(atm, 'c0') || isempty(atm.c0)
    atm.c0 = 20.05*sqrt(atm.Ts);
end
c0 = atm.c0;
T0 = atm.Ts;
a  = 11/6;                                   % von Karman exponent

fcol = f(:);
nf   = numel(fcol);

out.phi2 = zeros(nf, 3);  out.D = zeros(nf, 3);  out.omH = ones(nf, 3);
out.phi2max = inf(nf, 1); out.chi2weak = zeros(nf, 1); out.chi2 = zeros(nf, 1);
out.flim = NaN; out.Leff = 0;
sigma_dB = zeros(size(f));

if R1 <= 0 || hs <= hr
    return                                   % no slanted path (source not above receiver)
end

%% ---------------- (effective) propagation distance ----------------
cosTheta = (hs - hr)/R1;                     % theta from the vertical
if hs > atm.zi && opts.effectiveDistance
    fABL = max(atm.zi - hr, 0)/(hs - hr);    % fraction of the direct path below zi
    LABL = R1*fABL;
    Leff = LABL * (atm.zi/hs)^(-opts.q*(1 - atm.zi/hs));
    Leff = min(Leff, R1);
else
    Leff = R1;
end
out.Leff = Leff;
if Leff <= 0
    return
end
dh = Leff*cosTheta;                          % height span of the (effective) path

%% ---------------- integration grid along the path ----------------
% s = 1 - eta, from the receiver (s = 0) to the source end (s = 1). Dense
% near the receiver, where the length scales (~z) change fastest.
s = unique([linspace(0, 1, 400), logspace(-6, 0, 400)]).';
z = hr + dh*s;

T = get_abl_turbulence(z, atm);

Cg = sqrt(pi)*gamma(a-1)/gamma(a-1.5);

% strength-parameter densities along the path (integrands of LP24 Eqs. 35-36)
dT  = T.sigT2  .* T.LT  ./ T0^2;             % [nz x 1]
dvs = T.sigvs2 .* T.Lvs ./ c0^2;
dvb = T.sigvb2 * T.Lvb / c0^2 * ones(size(z));
IT  = trapz(s, dT);  Ivs = trapz(s, dvs);  Ivb = trapz(s, dvb);

%% ---------------- strength and wave parameters vs. frequency ----------------
k   = 2*pi*fcol/c0;
eta = 1 - s.';                               % [1 x nz], eta = 0 at the source end

phi2 = [Cg/2*k.^2*Leff*IT, 2*Cg*k.^2*Leff*Ivs, 2*Cg*k.^2*Leff*Ivb];   % [nf x 3]

% 1 - H(D(z), eta) with the LOCAL wave parameter D(z) = L/(k L_i(z)^2),
% inside the path integral as in LP24 Eq. (28); weighted with the strength
% density, so that <chi^2>_i = Phi_i^2/2 * omH_i exactly
Lscale = {T.LT, T.Lvs, T.Lvb*ones(size(z))};
dens   = {dT, dvs, dvb};
types  = {'T', 'v', 'v'};
omH = ones(nf, 3);  D = zeros(nf, 3);
for i = 1:3
    Dloc = Leff ./ (k * (Lscale{i}.').^2);                    % [nf x nz]
    c    = (eta.*(1 - eta)) .* Dloc;                          % argument of G
    G    = il_G(c, types{i});
    wgt  = dens{i}.';
    if trapz(s, dens{i}) > 0
        omH(:,i) = trapz(s, G .* wgt, 2) ./ trapz(s, dens{i});
    end
    D(:,i) = Leff ./ k * trapz(s, 1./max(Lscale{i}, realmin).^2);  % path-averaged D, output only
end

chi2_i   = phi2/2 .* omH;                    % per turbulence component, Eq. (28)
chi2weak = sum(chi2_i, 2);

phi2peak = max(phi2, [], 2);
beta     = phi2 ./ max(phi2peak, realmin);
phi2max  = 2*opts.b2 ./ max(sum(beta .* omH, 2), realmin);            % Eq. (30)

switch opts.saturationMode
    case 'total'        % DEFAULT, reproduces LKP25 Fig. 2 (see header)
        chi2 = chi2weak ./ (1 + phi2peak ./ phi2max);
    case 'component'    % LP24 Eq. (37) as printed
        chi2 = sum(chi2_i ./ (1 + phi2 ./ phi2max), 2);
    otherwise
        error('opts.saturationMode must be ''total'' or ''component''.')
end

if opts.saturation
    sigma_dB = reshape(20*log10(exp(1)) * sqrt(chi2), size(f));
else
    sigma_dB = reshape(20*log10(exp(1)) * sqrt(chi2weak), size(f));
end

idx = find(phi2peak >= phi2max, 1, 'first');
if ~isempty(idx), out.flim = fcol(idx); end

out.phi2 = phi2; out.D = D; out.omH = omH;
out.phi2max = phi2max; out.chi2weak = chi2weak; out.chi2 = chi2;

end

% =====================================================================
function G = il_G(c, type)
% G(c) = 1 - H for c = eta(1-eta)D, i.e. LP24 Eq. (22) (type 'T') or
% Eq. (25) (type 'v') written as Int_0^inf w(k)[1 - cos(c k)] dk.
% Interpolated (in log c) from a table computed once per session.
% c = Inf or NaN (at the ground for hr = 0) gives the limit G = 1.

persistent tab
if isempty(tab), tab = il_table(); end

if strcmp(type, 'T'), y = tab.GT; else, y = tab.Gv; end

G  = ones(size(c));
lc = log(c);
in = isfinite(lc) & lc >= tab.lc(1) & lc <= tab.lc(end);
G(in) = exp(interp1(tab.lc, log(y), lc(in), 'linear'));

lo = isfinite(lc) & lc < tab.lc(1);          % small c: power law of the first table interval
if any(lo(:))
    p = (log(y(2)) - log(y(1)))/(tab.lc(2) - tab.lc(1));
    G(lo) = y(1)*exp(p*(lc(lo) - tab.lc(1)));
end
G(c <= 0) = 0;
G = min(max(G, 0), 1);
end

% =====================================================================
function tab = il_table()
% Table of G(c) = Int_0^inf w(k)[1 - cos(c k)] dk on a log grid of c, for
%   w_T(k) = (a-1)/(1+k)^a,  w_v(k) = a(a-1) k/(1+k)^(a+1),  Int w dk = 1.
% The oscillatory part is integrated up to k0 = X/c on a fine log grid;
% beyond k0, 1 - cos is replaced by its mean 1 and the remaining weight is
% integrated analytically. G -> 1 for c -> Inf, G ~ c^(5/6) for c -> 0.

a  = 11/6;
X  = 1000;                                   % oscillation periods resolved explicitly ~ X/(2 pi)
tab.lc = linspace(log(1e-9), log(1e9), 361).';
GT = zeros(size(tab.lc)); Gv = GT;

for n = 1:numel(tab.lc)
    c  = exp(tab.lc(n));
    k0 = X/c;
    ds = 2*pi/(20*X);                        % 20 points per period at k0
    sg = (log(1e-12):ds:log(k0)).';
    kk = exp(sg);
    osc = 2*sin(c*kk/2).^2;                  % 1 - cos, numerically stable
    wT  = (a-1)./(1+kk).^a .* kk;            % includes dk = k ds
    wv  = a*(a-1)*kk./(1+kk).^(a+1) .* kk;
    GT(n) = trapz(sg, wT.*osc) + (1+k0)^(1-a);
    Gv(n) = trapz(sg, wv.*osc) + a*(1+k0)^(1-a) - (a-1)*(1+k0)^(-a);
end
tab.GT = max(GT, realmin);
tab.Gv = max(Gv, realmin);
end
