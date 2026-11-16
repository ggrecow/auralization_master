function map = get_doppler_mapping(t_e, tau)
% GET_DOPPLER_MAPPING  Mapping between emission and reception time of ONE sound path,
% used to apply the Doppler effect with a variable delay line.
%
%   map = get_doppler_mapping(t_e, tau)
%
% Sound emitted at the emission time t_e arrives at the receiver at the reception
% time t_r = t_e + tau(t_e), where tau is the propagation time of the eigenray.
% Reading the source signal at the emission time t_e(t_r) that belongs to each
% reception time gives the Doppler effect (frequency shift, and compression or
% stretching of the time axis) for any subsonic motion and any time-invariant
% medium, including refraction and wind (tau comes from the ray tracing).
%
% The Doppler factor is the ratio between reception and emission time increments:
%
%     D = dt_r/dt_e = 1 + dtau/dt_e     ( = 1 - M*cos(theta) for straight rays )
%
% so that a frequency f emitted by the source is received as f/D. D < 1 while the
% aircraft approaches, D > 1 while it moves away. D <= 0 (supersonic source
% towards the receiver) is not supported (error).
%
% PROCESSING
%   1) Smoothing of tau: the eigenray search has a finite accuracy (receiver
%      radius), which adds small random errors to tau. Since D is the derivative
%      of tau, these errors would become audible pitch and amplitude jitter. tau
%      is therefore smoothed with a local cubic least-squares fit (Savitzky-Golay
%      type) over a window of about 1 s. The fit is made on tau^2 and not on tau:
%      for a straight flight at constant speed in a homogeneous medium, tau^2 is
%      exactly a quadratic function of time (tau = R/c, R^2 = v^2*(t-t0)^2 + h^2),
%      so the smoothing does not change the physical propagation time, even close
%      to the receiver (a fit on tau would round off the minimum of tau at the
%      closest point). Accelerations, climbs and refraction add only slowly varying
%      terms (verified in <validate_doppler_delay_line>). The rms of the removed
%      part is returned in map.tau_smoothing_rms.
%   2) t_r = t_e + tau_smooth, checked to be strictly increasing (D > 0).
%   3) The inverse mapping t_e(t_r) is described by a cubic spline (smooth
%      Doppler factor, no steps at the time steps of the ray tracing).
%
% INPUTS
%   t_e : [K x 1] emission times of the ray-tracing time steps [s] (uniform dt)
%   tau : [K x 1] propagation time of the sound path at t_e [s]
%
% OUTPUT : struct
%   map.t_e               : [K x 1] emission times [s]
%   map.tau_raw           : [K x 1] propagation times, as given [s]
%   map.tau               : [K x 1] smoothed propagation times [s]
%   map.tau_smoothing_rms : rms of (tau_raw - tau) [s]
%   map.t_r               : [K x 1] reception times t_e + tau [s]
%   map.pp_te             : piecewise polynomial (spline) t_e(t_r), use with <ppval>
%   map.pp_dte            : piecewise polynomial of its derivative dt_e/dt_r = 1/D
%   map.D                 : [K x 1] Doppler factor dt_r/dt_e at t_e [-]
%
% See also: apply_doppler_delay_line, validate_doppler_delay_line
%
% Author: Gil Felix Greco, Braunschweig 16.11.2026
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to
%   develop this code under the author's direction. All codes were verified
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_e = t_e(:);
tau = tau(:);
K   = numel(t_e);

if numel(tau) ~= K
    error( 'auralization:dopplerSize', ...
        'get_doppler_mapping: <t_e> and <tau> must have the same number of elements.' );
end
if K < 4
    error( 'auralization:dopplerSize', ...
        'get_doppler_mapping: at least 4 time steps are needed (got %d).', K );
end

%% 1) smoothing of the propagation time (local cubic fit of tau^2 over approx. 1 s)

dt = median( diff(t_e) );
t_window = 1.0;                                 % [s] smoothing window
nHalf = max( 2, round( 0.5*t_window/dt ) );     % half window, in time steps
nWin  = min( K, 2*nHalf + 1 );
polyOrder = 3;

tau2 = tau.^2;
tau2_s = zeros(K,1);
for i = 1:K
    i1 = max( 1, i - nHalf );
    i2 = min( K, i1 + nWin - 1 );
    i1 = max( 1, i2 - nWin + 1 );               % window shifted at the edges, same length
    idx = (i1:i2).';
    tt = ( t_e(idx) - t_e(i) ) / dt;            % local time, scaled for conditioning
    A  = tt .^ (0:polyOrder);
    c  = A \ tau2(idx);
    tau2_s(i) = c(1);
end
tau_s = sqrt( max( tau2_s, 0 ) );

%% 2) reception times

t_r = t_e + tau_s;

if any( diff(t_r) <= 0 )
    error( 'auralization:dopplerNotMonotonic', ...
        ['get_doppler_mapping: the reception times t_e + tau do not increase monotonically ' ...
         '(Doppler factor <= 0). This happens for a source moving towards the receiver faster ' ...
         'than sound, or for erroneous propagation times. The delay-line Doppler method ' ...
         '(doppler_method = delay_line) cannot be used for this case.'] );
end

%% 3) inverse mapping t_e(t_r) and its derivative

pp_te = spline( t_r, t_e );
pp_dte = il_ppder( pp_te );

D = 1 ./ ppval( pp_dte, t_r );

%% output

map.t_e = t_e;
map.tau_raw = tau;
map.tau = tau_s;
map.tau_smoothing_rms = sqrt( mean( (tau - tau_s).^2 ) );
map.t_r = t_r;
map.pp_te = pp_te;
map.pp_dte = pp_dte;
map.D = D;

end

%% local function : derivative of a piecewise polynomial (no toolbox needed)
function ppd = il_ppder(pp)
[breaks, coefs, nPieces, order, dim] = unmkpp(pp);
if order <= 1
    ppd = mkpp( breaks, zeros(nPieces*dim, 1), dim );
    return
end
powers = (order-1):-1:1;
dcoefs = coefs(:, 1:end-1) .* powers;
ppd = mkpp( breaks, dcoefs, dim );
end
