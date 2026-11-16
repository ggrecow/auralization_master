function [y, gain] = apply_doppler_delay_line(x, fs, map, t_out, convective_amplification)
% APPLY_DOPPLER_DELAY_LINE  Doppler effect (and convective amplification) of ONE sound
% path using a variable delay line.
%
%   [y, gain] = apply_doppler_delay_line(x, fs, map, t_out, convective_amplification)
%
% The source signal x(t_e), synthesized in emission time from source data WITHOUT
% Doppler effect, is read at the emission time that belongs to each reception time:
%
%     y(t_r) = G(t_r) * x( t_e(t_r) ),     t_r = t_e + tau(t_e)
%
% where tau is the propagation time of the sound path (see <get_doppler_mapping>).
% This produces the Doppler frequency shift f/D and the compression (approach) or
% stretching (departure) of the time axis, with D = dt_r/dt_e the Doppler factor.
% The propagation filters (spreading, absorption, ground reflection, HRTFs) are
% applied AFTER this step, in reception time, i.e. at the received frequencies.
%
% Fractional-delay interpolation: Blackman-windowed sinc, 32 taps. Where the
% signal is read faster than it was sampled (D < 1, approach), the cut-off of the
% kernel is lowered to D*fs/2 to avoid aliasing (band-limited resampling).
% Accuracy (see <validate_doppler_delay_line>): level error < 0.01 dB up to 0.35*fs.
%
% Convective amplification (optional): a source moving towards the receiver is
% louder than the same source at rest, by the factor
%
%     G = D^(-2)  on the sound pressure,  i.e.  dL = -40*log10(D)  [dB]
%
% with D = 1 - M*cos(theta) for straight rays (classical result for a moving point
% source, e.g. Dowling & Ffowcs Williams 1983, Ch. 9; the same 40*log10 term is used
% in aircraft noise prediction to project static source levels to flight). The delay
% line itself keeps the amplitude of the signal; G is only applied if
% <convective_amplification> = 1. Use 0 if the source data already contain it.
% D is taken from the propagation times, so it is consistent with the Doppler shift
% of each sound path.
%
% INPUTS
%   x     : [N x 1] source signal sampled at t_e = (n-1)/fs, n = 1..N [Pa]
%   fs    : sampling frequency [Hz]
%   map   : struct from <get_doppler_mapping> (one sound path)
%   t_out : [M x 1] reception times of the output samples [s]. The corresponding
%           emission times must lie within the signal x (checked).
%   convective_amplification : 0 or 1 (DEFAULT 1)
%
% OUTPUTS
%   y    : [M x 1] signal at the receiver side of the path, before the
%          propagation filters [Pa]
%   gain : [M x 1] convective amplification factor applied (ones if off)
%
% See also: get_doppler_mapping, validate_doppler_delay_line, apply_propagation
%
% Author: Gil Felix Greco, Braunschweig 16.11.2026
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to
%   develop this code under the author's direction. All codes were verified
%   by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5 || isempty(convective_amplification)
    convective_amplification = 1;
end

x = x(:);
t_out = t_out(:);
N = numel(x);

%% emission time and Doppler factor at each output sample

te_out = ppval( map.pp_te, t_out );         % emission time heard at t_out [s]
D_out  = 1 ./ ppval( map.pp_dte, t_out );   % Doppler factor dt_r/dt_e [-]

tol = 1e-9;
if min(te_out) < -tol || max(te_out) > (N-1)/fs + tol
    error( 'auralization:dopplerRange', ...
        ['apply_doppler_delay_line: the requested reception times need emission times between ' ...
         '%.4f s and %.4f s, but the source signal covers 0 to %.4f s.'], ...
         min(te_out), max(te_out), (N-1)/fs );
end

%% band-limited fractional-delay interpolation (Blackman-windowed sinc)

half = 16;                                  % 2*half = 32 taps
pos  = te_out*fs + 1;                       % fractional read index into x (1-based)
n0   = floor(pos);
d    = pos - n0;                            % fractional part, 0 <= d < 1
fc   = min( 1, D_out );                     % normalized cut-off (anti-aliasing if D < 1)

xp = [ zeros(half,1); x; zeros(half+1,1) ]; % zero padding: index n in x -> n + half in xp
y  = zeros( size(pos) );

for k = -half+1:half
    a = k - d;                              % distance to the tap [samples]
    w = 0.42 + 0.5*cos( pi*a/half ) + 0.08*cos( 2*pi*a/half );   % Blackman window
    s = fc .* il_sinc( fc .* a );
    y = y + xp( n0 + k + half ) .* s .* w;
end

%% convective amplification

if convective_amplification == 1
    gain = D_out .^ (-2);
    y = y .* gain;
else
    gain = ones( size(y) );
end

end

%% local function : normalized sinc, sin(pi*a)/(pi*a), no toolbox needed
function s = il_sinc(a)
s = ones( size(a) );
nz = a ~= 0;
s(nz) = sin( pi*a(nz) ) ./ ( pi*a(nz) );
end
