function validate_time_step( PATH_source )
% function validate_time_step( PATH_source )
%
% Verification of <get_time_step> (time step of the input data estimated from all time steps,
% rounded to an integer number of samples, consistency check).
%
%   1) synthetic time vectors: exact and rounded (0.1 s) time values of a 0.125 s time step,
%      missing and repeated time steps, a varying time step, random jitter. Checks the time step
%      found and whether the consistency check passes (no warning) or fails (warning)
%   2) example PANAM data (input_data/auralization_input.dat, or <PATH_source>): the source time is
%      rounded to 0.1 s (differences of 0.1 and 0.2 s), while the time step is 0.125 s. Checks
%      dt = 0.125 s (6000 samples at 48 kHz) for the full record and for parts of it (as after
%      trimming), and compares with two estimates from other quantities of the data
%      (distance per step / vel, and (retarded time - dist/c) per step)
%
% INPUT (optional):
%   PATH_source : string
%   path of a PANAM input file (DEFAULT: input_data/auralization_input.dat of this repository)
%
% Gil Felix Greco, Braunschweig 18.11.2026
%
%   AI disclosure: Generative AI assistance (Anthropic Claude) was used to develop this
%   code under the author's direction. All codes were verified by the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1 || isempty( PATH_source )
    PATH_source = fullfile( fileparts( fileparts( mfilename('fullpath') ) ), 'input_data', 'auralization_input.dat' );
end

fs  = 48000;
dt0 = 0.125;
nFail = 0;

fprintf('\n*--------------------------------------------------------------------------*\n');
fprintf('validate_time_step\n');

%% 1) synthetic time vectors

q = @(x) round( x*10 )/10;   % rounding to 0.1 s
k = ( 0:160 ).';
rng( 1 );

% name, time values, expected number of samples (NaN: not checked), expected consistent
cases = {
    'exact, uniform',                       k*dt0 + 3,                                   6000,  true;
    'rounded to 0.1 s',                     q( k*dt0 + dt0 ),                            6000,  true;
    'rounded, 1 missing step (middle)',     q( [k(1:80); k(82:end)]*dt0 + dt0 ),         NaN,   false;
    'rounded, 1 missing step (start)',      q( [k(1:5); k(7:end)]*dt0 + dt0 ),           NaN,   false;
    'rounded, repeated step (+1 ms)',       sort( [q( k*dt0 + dt0 ); q( 10*dt0 + dt0 ) + 1e-3] ), NaN, false;
    'varying step 0.10 -> 0.15 s',          cumsum( [0; linspace( 0.10, 0.15, 160 ).'] ), NaN,  false;
    'random jitter, 10 ms rms',             k*dt0 + 0.01*randn( 161, 1 ),                NaN,   true;
    'dt = 0.5 s, rounded to 0.1 s',         q( k*0.5 + 0.05 ),                           24000, true;
    };

fprintf('\n1) synthetic time vectors (fs = %d Hz)\n', fs);
for c = 1:size( cases, 1 )
    lastwarn('');
    s = warning( 'off', 'auralization:timeStepNotUniform' );
    [dt, info] = get_time_step( cases{c,2}, fs, 'validate_time_step' );
    warning( s );
    ok = ( info.consistent == cases{c,4} ) && ( isnan( cases{c,3} ) || info.nSamples == cases{c,3} );
    nFail = nFail + ~ok;
    fprintf('  %-36s dt = %.6f s (%5d samples), max. deviation %.4f s (limit %.4f s), consistent = %d  -> %s\n', ...
        cases{c,1}, dt, info.nSamples, info.max_deviation, dt/2, info.consistent, il_ok( ok ) );
end

% 44.1 kHz: 0.125 s = 5512.5 samples -> 5510 (multiple of 10 samples)
[~, info] = get_time_step( q( k*dt0 + dt0 ), 44100, 'validate_time_step' );
ok = info.nSamples == 5510 && info.consistent;
nFail = nFail + ~ok;
fprintf('  %-36s dt = %.6f s (%5d samples)  -> %s\n', 'rounded to 0.1 s, fs = 44.1 kHz', info.dt, info.nSamples, il_ok( ok ) );

% errors
try
    get_time_step( [0; 0.1; 0.1; 0.2], fs, 'validate_time_step' ); ok = false;
catch
    ok = true;
end
nFail = nFail + ~ok;
fprintf('  %-36s -> %s\n', 'repeated time value: error', il_ok( ok ) );

%% 2) example PANAM data

fprintf('\n2) example PANAM data: %s\n', PATH_source);
if ~isfile( PATH_source )
    fprintf('  file not found, skipped\n');
else
    data = PANAM2struct( PATH_source );
    data = data(:,1);
    K = numel( data );
    t = zeros( K, 1 ); te = t; pos = zeros( K, 3 ); vel = t;
    for i = 1:K
        t(i)     = data{i}.source_time;
        te(i)    = data{i}.retarded_time - data{i}.dist/data{i}.sound_speed;
        pos(i,:) = [ data{i}.xpos data{i}.ypos data{i}.alt ];
        vel(i)   = data{i}.vel;
    end

    [dt, info] = get_time_step( t, fs, 'validate_time_step' );
    ok = info.nSamples == 6000 && info.consistent;
    nFail = nFail + ~ok;
    fprintf('  full record (%d time steps): dt = %.6f s (%d samples), first difference %.4f s, differences %.4f ... %.4f s, max. deviation %.4f s  -> %s\n', ...
        K, dt, info.nSamples, info.dt_first, info.diff_min, info.diff_max, info.max_deviation, il_ok( ok ) );
    fprintf('  for comparison: (last - first)/(K - 1) = %.6f s; distance per step / vel = %.5f s; (retarded time - dist/c) per step = %.5f s\n', ...
        ( t(end) - t(1) )/(K-1), sum( sqrt( sum( diff(pos).^2, 2 ) ) )/(K-1)/mean(vel), ( te(end) - te(1) )/(K-1) );

    % parts of the record (as after trimming): all lengths >= 12 time steps, every 5th start
    nParts = 0; nBad = 0;
    for W = 12:K
        for i0 = 1:5:K-W+1
            [~, info] = get_time_step( t(i0:i0+W-1), fs, 'validate_time_step' );
            nParts = nParts + 1;
            nBad = nBad + ~( info.nSamples == 6000 && info.consistent );
        end
    end
    ok = nBad == 0;
    nFail = nFail + ~ok;
    fprintf('  %d parts of the record (12 to %d time steps): %d with dt ~= 6000 samples or failed check  -> %s\n', ...
        nParts, K, nBad, il_ok( ok ) );
end

%% summary

if nFail == 0
    fprintf('\nvalidate_time_step: all checks OK\n');
else
    fprintf('\nvalidate_time_step: %d check(s) FAILED\n', nFail);
end
fprintf('*--------------------------------------------------------------------------*\n');

end

function s = il_ok( ok )
if ok
    s = 'OK';
else
    s = 'FAILED';
end
end
