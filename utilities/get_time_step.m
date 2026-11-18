function [dt, info] = get_time_step( time, fs, caller, nMultiple )
% function [dt, info] = get_time_step( time, fs, caller, nMultiple )
%
% Time step (dt) of the sound emission predictions (e.g. PANAM), estimated from the time values
% of ALL time steps, and check that the time values fit a uniform time grid.
%
% The auralization needs a constant time step: each time step of the input data is rendered as
% one block of fs*dt samples (broadband synthesis, tonal envelopes, propagation filters, coherence
% loss, delay lines). The time values in the input data can be rounded: PANAM writes the source
% time with a resolution of 0.1 s, while its time step can be e.g. 0.125 s. The differences
% between consecutive time values are then 0.1 or 0.2 s, and the first difference is not the
% time step. This function therefore does not use the first difference:
%
%   1) dt is the step of the uniform grid t0 + (k-1)*dt that is closest to all time values
%      (smallest maximum deviation, i.e. minimax line fit). Rounding errors do not accumulate:
%      for time values rounded to a resolution q, dt is off by at most q/(K-1), K = number of
%      time values, and in practice much less, as the rounding errors of all values constrain
%      the fit (example PANAM data, q = 0.1 s: dt = 0.125 s exactly for the full record and for
%      every part of it with 12 or more time steps).
%   2) dt is rounded so that the block length fs*dt is an integer MULTIPLE OF <nMultiple> samples
%      (DEFAULT 10), dt = nMultiple*round(fs*dt/nMultiple)/fs. The broadband synthesis
%      (<broadbandSynthesis_smooth>) splits each block into 5 hops (HopSizeFactor = 5) and needs
%      an even block length for its one-sided spectrum, so the block length must be a multiple
%      of 10 samples. An integer number of samples alone is not enough: e.g. 4802 samples
%      (dt = 0.100042 s at 48 kHz) gives a hop size of 960.4 samples and the synthesis fails.
%      Examples: 0.125 s -> 6000 samples at 48 kHz; 0.1000337 s -> 4800 samples (dt = 0.1 s,
%      0.03 % faster); at 44.1 kHz, 0.125 s = 5512.5 -> 5510 samples, i.e. dt = 0.124943 s
%      (0.05 % faster).
%   3) Consistency check: every time value must lie within dt/2 of the uniform grid (offset t0
%      fitted, not taken from the first value). Then each time step of the input data falls into
%      its own block (none skipped or doubled). If the check fails (e.g. missing time steps,
%      repeated time steps not removed by <prepare_input_SQ>, or a time step that really varies),
%      a WARNING names the time value with the largest deviation. The run continues, with each
%      time step rendered as one block of dt.
%      Limitation: a single missing time step in an otherwise EXACT time column cannot be told
%      apart from a slightly longer dt (the deviations then stay just below dt/2); with rounded
%      time values (as in PANAM), it exceeds the limit and is detected.
%
% The input data are not changed: spectra, positions, angles and the order of the time steps are
% kept. Only the time assigned to each time step changes, from the (rounded) value in the data to
% (k-1)*dt (see <getAuralizationTime>).
%
% INPUTS:
%   time : vector [K x 1] or [1 x K]
%   time values of the time steps (e.g. <source_time> of the input data), in seconds, increasing
%
%   fs : scalar (optional)
%   sampling frequency used for the auralization [Hz]. If empty or not given, <sampling_freq>
%   from the .ini file is used (DEFAULT 48000 Hz, as in <MASTER_auralization_EngineAirframe>)
%
%   caller : string (optional)
%   name of the calling function, used in the warning
%
%   nMultiple : positive integer (optional)
%   the block length fs*dt is rounded to a multiple of nMultiple samples. DEFAULT: 10, required
%   by <broadbandSynthesis_smooth> (HopSizeFactor = 5, even block length). Use 1 only if the
%   broadband synthesis is not used.
%
% OUTPUTS:
%   dt : scalar
%   time step [s], an integer multiple of nMultiple samples at fs
%
%   info : struct
%   info.dt            : dt [s] (same as output dt)
%   info.dt_fit        : dt before rounding to an integer number of samples [s]
%   info.nSamples      : block length fs*dt [samples], a multiple of info.nMultiple
%   info.nMultiple     : nMultiple used [samples]
%   info.fs            : sampling frequency used [Hz]
%   info.K             : number of time values
%   info.dt_first      : first difference of the time values, time(2)-time(1) [s] (for comparison)
%   info.diff_min      : smallest difference between consecutive time values [s]
%   info.diff_max      : largest difference between consecutive time values [s]
%   info.t0            : fitted time of the first grid point [s]
%   info.max_deviation : largest deviation of the time values from the grid [s]
%   info.idx_worst     : index of the time value with the largest deviation
%   info.consistent    : true if info.max_deviation < dt/2
%
% Gil Felix Greco, Braunschweig 18.11.2026
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to develop this
%   code under the author's direction. All codes were verified by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global input_file

if nargin < 3 || isempty( caller )
    caller = 'get_time_step';
end

if nargin < 4 || isempty( nMultiple )
    nMultiple = 10; % block length must be a multiple of 2*HopSizeFactor (see <broadbandSynthesis_smooth>)
end

if ~isscalar( nMultiple ) || nMultiple < 1 || nMultiple ~= round( nMultiple )
    error( 'auralization:timeStep', '%s: nMultiple must be a positive integer.', caller );
end

if nargin < 2 || isempty( fs )
    if isstruct( input_file ) && isfield( input_file, 'sampling_freq' )
        fs = str2double( input_file.sampling_freq );
    else
        fs = 48000; % default value (as in <MASTER_auralization_EngineAirframe>)
    end
end

if ~isscalar( fs ) || ~isfinite( fs ) || fs <= 0
    error( 'auralization:invalidSamplingFrequency', ...
        '<sampling_freq> in the .ini file must be a positive number.' );
end

time = double( time(:) );
K = numel( time );

if K < 2 || any( ~isfinite( time ) )
    error( 'auralization:timeStep', ...
        '%s: at least 2 finite time values are needed to determine the time step (got %d).', caller, K );
end

d = diff( time );
if any( d <= 0 )
    idx = find( d <= 0, 1 );
    error( 'auralization:timeStep', ...
        ['%s: the time values of the input data must increase (time value %d: %.4f s, ' ...
         'time value %d: %.4f s).'], caller, idx, time(idx), idx+1, time(idx+1) );
end

k = ( 0:K-1 ).';

%% 1) step of the uniform grid closest to all time values (minimax line fit)
% For a given step b, the best offset puts the grid midway between the largest and smallest
% residual, so the maximum deviation is f(b) = ( max(r) - min(r) )/2, r = time - k*b.
% f(b) is convex and piecewise linear; its minimum lies between the smallest and largest
% difference of consecutive time values. It is found by bisection on the sign of the slope of f.

lo = min( d );
hi = max( d );
for iter = 1:200
    if hi - lo <= 4*eps( hi )
        break
    end
    b = ( lo + hi )/2;
    r = time - k*b;
    [~, iMax] = max( r );
    [~, iMin] = min( r );
    slope = ( k(iMin) - k(iMax) )/2; % slope of f(b)
    if slope > 0
        hi = b;
    elseif slope < 0
        lo = b;
    else
        lo = b; hi = b;
    end
end
dt_fit = ( lo + hi )/2;

%% 2) integer multiple of nMultiple samples

% (fs*dt/nMultiple is first rounded to 1e-6, so that a value of exactly x.5 is always rounded up,
% independent of floating-point errors of the fit)
nSamples = nMultiple*round( round( fs*dt_fit/nMultiple*1e6 )/1e6 );
if nSamples < nMultiple
    error( 'auralization:timeStep', ...
        '%s: the time step dt = %g s is shorter than %d samples at fs = %g Hz.', caller, dt_fit, nMultiple, fs );
end
dt = nSamples/fs;

%% 3) consistency check (with the dt actually used)

r  = time - k*dt;
t0 = ( max(r) + min(r) )/2;
[max_deviation, idx_worst] = max( abs( r - t0 ) );
consistent = max_deviation < dt/2;

if ~consistent
    warning( 'auralization:timeStepNotUniform', ...
        ['%s: the time values of the input data do not fit a uniform time step. ' ...
         'dt = %.6f s (from all %d time values). Time value %d (t = %.4f s) deviates by %.4f s from ' ...
         'the uniform grid (limit dt/2 = %.4f s). Differences between consecutive time values: ' ...
         '%.4f ... %.4f s. Possible causes: missing or repeated time steps, or a time step that ' ...
         'really varies. Each time step is still rendered as one block of dt: check the input data.'], ...
        caller, dt, K, idx_worst, time(idx_worst), max_deviation, dt/2, min(d), max(d) );
end

%% output

info.dt            = dt;
info.dt_fit        = dt_fit;
info.nSamples      = nSamples;
info.nMultiple     = nMultiple;
info.fs            = fs;
info.K             = K;
info.dt_first      = d(1);
info.diff_min      = min( d );
info.diff_max      = max( d );
info.t0            = t0;
info.max_deviation = max_deviation;
info.idx_worst     = idx_worst;
info.consistent    = consistent;

end
